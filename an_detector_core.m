function out = an_detector_core(mode, varargin)
%AN_DETECTOR_CORE Lean detector-map utility for the PILATUS LUT package.
%
% This package-specific core intentionally contains only q-phi binning.
% Detector-pose calibration is performed by fit_pilatus100k_lut_robust.m.
% It does NOT call or require estSharedRingOrigin.m, cif2powder_CM.m, or
% LDSDw_all.m.
%
% Public modes
% ------------
%   cfg = an_detector_core('defaults')
%       Return a configuration structure containing cfg.tripol.
%
%   opts = an_detector_core('tripol_defaults')
%       Return the q-phi options directly.
%
%   TP = an_detector_core('tripol', image, geom, opts)
%       Bin an image on the native detector grid into fixed and/or adaptive
%       q-phi representations. opts may be omitted or partially specified.
%
% Required geometry fields
% ------------------------
%   geom.qMap       momentum transfer at every detector pixel
%   geom.phiMap     azimuth at every detector pixel, in radians
%   geom.totalCorr  required only when opts.applyTotalCorr is true
%
% Output
% ------
%   TP.mode1        fixed q-phi grid
%   TP.mode2        adaptive phi bins for each q ring
%   TP.fitReady     finite weighted values from adaptive mode
%   TP.meta         options and bookkeeping
%
% This file replaces the older consolidated detector core, whose
% 'analyze_lab6' mode depended on estSharedRingOrigin.m. That mode is not
% part of this LUT package.
%
% Adi Natan, 2026.

if nargin < 1 || isempty(mode)
    error(['Specify a mode: an_detector_core(''tripol'', image, geom, opts), ' ...
        'an_detector_core(''defaults''), or an_detector_core(''tripol_defaults'').']);
end

mode = lower(char(mode));
switch mode
    case 'defaults'
        out = struct('tripol', local_default_opts());

    case {'tripol_defaults','tripol-defaults','tripoldefaults'}
        out = local_default_opts();

    case 'tripol'
        if numel(varargin) < 2
            error('Use an_detector_core(''tripol'', image, geom, opts).');
        end
        image = varargin{1};
        geom = varargin{2};
        if numel(varargin) >= 3 && ~isempty(varargin{3})
            opts = varargin{3};
        else
            opts = struct();
        end
        out = local_tripol_qphi(image, geom, opts);

    case 'info'
        out = struct( ...
            'name', 'an_detector_core', ...
            'purpose', 'package-specific fast q-phi binning', ...
            'modes', {{'defaults','tripol_defaults','tripol'}}, ...
            'externalDependencies', {{}}, ...
            'supportsAnalyzeLab6', false);

    case 'analyze_lab6'
        error(['The package-specific an_detector_core no longer implements ' ...
            '''analyze_lab6''. Use fit_pilatus100k_lut_robust.m for detector ' ...
            'pose calibration. estSharedRingOrigin.m is not required.']);

    otherwise
        error('Unknown mode ''%s''. Use defaults, tripol_defaults, tripol, or info.', mode);
end
end

function out = local_tripol_qphi(im2, geom, optsIn)
%LOCAL_TRIPOL_QPHI Fast q-phi binning on the native detector grid.
%
%   out = local_tripol_qphi(im2, geom, opts)
%
% The output fields match the package q-phi interface:
%   out.mode1     fixed q-phi grid
%   out.mode2     adaptive phi grid for each q ring
%   out.fitReady  finite, weighted adaptive-bin values
%   out.meta
%
% Speed comes from processing only occupied bins.  It avoids the original
% nested loop over every possible q-by-phi bin while preserving trimmean,
% MAD-SEM, variance propagation, SEM floors, and the continuous requested
% phi window.

if nargin < 3 || isempty(optsIn)
    optsIn = struct();
end
opts = local_merge_struct(local_default_opts(), optsIn);

if ~isfield(geom,'qMap') || ~isfield(geom,'phiMap')
    error('geom must contain qMap and phiMap.');
end
if ~isequal(size(im2), size(geom.qMap), size(geom.phiMap))
    error('im2, geom.qMap, and geom.phiMap must have identical sizes.');
end

qEdges = opts.qEdges(:).';
if numel(qEdges) < 2 || any(diff(qEdges) <= 0)
    error('opts.qEdges must be monotonically increasing.');
end
qCenters = qEdges(1:end-1) + 0.5*diff(qEdges);
nQ = numel(qCenters);

phiEdgesFixed = opts.fixedPhiEdges(:).';
if numel(phiEdgesFixed) < 2 || any(diff(phiEdgesFixed) <= 0)
    error('opts.fixedPhiEdges must be monotonically increasing.');
end
phiLo = phiEdgesFixed(1);
phiHi = phiEdgesFixed(end);
phiWidth = phiHi - phiLo;
if phiWidth <= 0 || phiWidth > 2*pi + 100*eps
    error('The fixed phi window must be positive and no wider than 2*pi.');
end
phiCtr = 0.5*(phiLo + phiHi);

if opts.applyTotalCorr
    if ~isfield(geom,'totalCorr')
        error('geom.totalCorr is required when applyTotalCorr is true.');
    end
    corrMap = geom.totalCorr;
else
    corrMap = ones(size(im2));
end
Icorr = im2 .* corrMap;
phiMapUse = local_wrap_to_center(geom.phiMap - deg2rad(opts.phi0_deg), phiCtr);

mask = isfinite(im2) & isfinite(Icorr) & isfinite(geom.qMap) & isfinite(phiMapUse);
if ~isempty(opts.mask)
    if ~isequal(size(opts.mask), size(im2))
        error('opts.mask must have the same size as im2.');
    end
    mask = mask & opts.mask;
end
if ~isempty(opts.minIcorr)
    mask = mask & Icorr >= opts.minIcorr;
end
if ~isempty(opts.maxIcorr)
    mask = mask & Icorr <= opts.maxIcorr;
end
mask = mask & geom.qMap >= qEdges(1) & geom.qMap <= qEdges(end) & ...
    phiMapUse >= phiLo & phiMapUse <= phiHi;

qVals = geom.qMap(mask);
phiVals = phiMapUse(mask);
Ivals = Icorr(mask);

hasVarMap = ~isempty(opts.varMap);
if hasVarMap
    if ~isequal(size(opts.varMap), size(im2))
        error('opts.varMap must have the same size as im2.');
    end
    if opts.applyTotalCorr
        Vvals = opts.varMap(mask) .* corrMap(mask).^2;
    else
        Vvals = opts.varMap(mask);
    end
else
    Vvals = [];
end

qbin = discretize(qVals, qEdges);
valid = isfinite(qbin) & qbin >= 1 & qbin <= nQ & isfinite(Ivals);
qbin = qbin(valid);
phiVals = phiVals(valid);
Ivals = Ivals(valid);
if hasVarMap
    Vvals = Vvals(valid);
end

% Sorted q-bin ranges are reused by adaptive mode.
[qbinSorted, ordQ] = sort(qbin);
nValid = numel(qbinSorted);
firstQ = zeros(nQ,1);
lastQ = zeros(nQ,1);
if nValid > 0
    firstQ = accumarray(qbinSorted, (1:nValid).', [nQ 1], @min, 0);
    lastQ = accumarray(qbinSorted, (1:nValid).', [nQ 1], @max, 0);
end

if opts.computeMode2
    mode2 = local_build_mode2(qCenters, qEdges, phiVals, Ivals, Vvals, ...
        hasVarMap, ordQ, firstQ, lastQ, phiLo, phiHi, opts);
else
    mode2 = local_empty_mode2(qCenters, qEdges);
end

if opts.computeMode1
    mode1 = local_build_mode1(qCenters, qEdges, phiEdgesFixed, qbin, ...
        phiVals, Ivals, Vvals, hasVarMap, opts);
else
    mode1 = local_empty_mode1(qCenters, qEdges, phiEdgesFixed);
end

fitReady = local_make_fit_ready(mode2, opts);

out = struct();
out.mode2 = mode2;
out.mode1 = mode1;
out.fitReady = fitReady;
out.meta = struct();
out.meta.optsUsed = opts;
out.meta.validPixelsUsed = numel(Ivals);
out.meta.applyTotalCorr = opts.applyTotalCorr;
out.meta.phi0_deg = opts.phi0_deg;
out.meta.hasVarMap = hasVarMap;
out.meta.phiEdgesUsed = phiEdgesFixed(:);
out.meta.phiWindow_deg = rad2deg([phiLo phiHi]);
out.meta.phiCenter_deg = rad2deg(phiCtr);
out.meta.mode2UsesFixedPhiWindow = true;
out.meta.engine = 'an_detector_core:tripol-fast';
end

% ======================================================================
function opts = local_default_opts()
opts = struct();
opts.qEdges = 0.3:0.0025:9.1;
opts.fixedPhiEdges = linspace(-pi, pi, 361);
opts.applyTotalCorr = true;
opts.mask = [];
opts.varMap = [];
opts.phi0_deg = 0;
opts.targetPixPerBin = 80;
opts.minPhiBins = 16;
opts.maxPhiBins = 256;
opts.minPixPerBin = 6;
opts.minPixPerRing = 60;
opts.forceEvenPhiBins = true;
opts.binMean = 'trimmean';
opts.trimPct = 37;
opts.errModel = 'madsem';
opts.semFloorFrac = 0.25;
opts.semFloorAbs = 0;
opts.minIcorr = [];
opts.maxIcorr = [];
opts.computeMode1 = true;
opts.computeMode2 = true;
end

% ======================================================================
function mode2 = local_build_mode2(qCenters, qEdges, phiVals, Ivals, Vvals, ...
    hasVarMap, ordQ, firstQ, lastQ, phiLo, phiHi, opts)

nQ = numel(qCenters);
mode2 = local_empty_mode2(qCenters, qEdges);

for iq = 1:nQ
    if firstQ(iq) == 0
        continue
    end
    idx = ordQ(firstQ(iq):lastQ(iq));
    nq = numel(idx);
    mode2.Nq(iq) = nq;
    if nq < opts.minPixPerRing
        continue
    end

    Iqv = Ivals(idx);
    if hasVarMap
        Vqv = Vvals(idx);
    else
        Vqv = [];
    end
    [mode2.Iq(iq), mode2.IqSem(iq), mode2.IqStd(iq)] = ...
        local_bin_stats(Iqv, Vqv, opts);

    Nphi = local_choose_nphi(nq, opts);
    mode2.Nphi(iq) = Nphi;
    phiEdges = linspace(phiLo, phiHi, Nphi+1);
    phiCenters = phiEdges(1:end-1) + 0.5*diff(phiEdges);

    % Uniform-bin arithmetic is faster than histcounts for every q ring.
    rel = (phiVals(idx) - phiLo) / (phiHi - phiLo);
    pbin = floor(rel*Nphi) + 1;
    pbin(rel >= 1) = Nphi;
    pbin = min(max(pbin,1),Nphi);

    [mu, sem, stdv, np] = local_grouped_stats(Iqv, Vqv, pbin, Nphi, opts);
    sem = local_apply_sem_floor(sem, np, opts);
    wPhi = nan(size(sem));
    goodW = isfinite(sem) & sem > 0 & isfinite(mu);
    wPhi(goodW) = 1 ./ sem(goodW).^2;

    mode2.phiEdges{iq} = phiEdges(:);
    mode2.phiCenters{iq} = phiCenters(:);
    mode2.Imean{iq} = mu;
    mode2.Isem{iq} = sem;
    mode2.Istd{iq} = stdv;
    mode2.Npix{iq} = np;
    mode2.wPhi{iq} = wPhi;
end
end

% ======================================================================
function mode1 = local_build_mode1(qCenters, qEdges, phiEdges, qbin, ...
    phiVals, Ivals, Vvals, hasVarMap, opts)

nQ = numel(qCenters);
phiCenters = phiEdges(1:end-1) + 0.5*diff(phiEdges);
nPhi = numel(phiCenters);
mode1 = local_empty_mode1(qCenters, qEdges, phiEdges);

pbin = discretize(phiVals, phiEdges);
valid = isfinite(pbin) & pbin >= 1 & pbin <= nPhi;
linBin = qbin(valid) + (pbin(valid)-1)*nQ;
if hasVarMap
    Vuse = Vvals(valid);
else
    Vuse = [];
end
[mu, sem, stdv, np] = local_grouped_stats(Ivals(valid), Vuse, linBin, nQ*nPhi, opts);

mode1.Iqp = reshape(mu, nQ, nPhi);
mode1.Semqp = reshape(sem, nQ, nPhi);
mode1.Stdqp = reshape(stdv, nQ, nPhi);
mode1.Nqp = reshape(np, nQ, nPhi);

for iq = 1:nQ
    mode1.Semqp(iq,:) = local_apply_sem_floor(mode1.Semqp(iq,:).', ...
        mode1.Nqp(iq,:).', opts).';
end
mode1.Wqp = nan(nQ,nPhi);
goodW = isfinite(mode1.Semqp) & mode1.Semqp > 0 & isfinite(mode1.Iqp);
mode1.Wqp(goodW) = 1 ./ mode1.Semqp(goodW).^2;
end

% ======================================================================
function [mu, sem, stdv, counts] = local_grouped_stats(vals, vars, groups, nGroups, opts)
vals = vals(:);
groups = groups(:);
hasVars = ~isempty(vars);
if hasVars
    vars = vars(:);
end

mu = nan(nGroups,1);
sem = nan(nGroups,1);
stdv = nan(nGroups,1);
counts = accumarray(groups, ones(size(groups)), [nGroups 1], @sum, 0);
if isempty(groups)
    return
end

[gSorted, ord] = sort(groups);
starts = [1; find(diff(gSorted) ~= 0) + 1];
stops = [starts(2:end)-1; numel(gSorted)];
groupIDs = gSorted(starts);

nOcc = stops-starts+1;
validOcc = nOcc >= opts.minPixPerBin;

% Single-pixel groups are common on fine fixed grids.  Handle all of them
% at once instead of calling the robust-statistics helper thousands of times.
single = validOcc & nOcc == 1;
if any(single)
    gid = groupIDs(single);
    id = ord(starts(single));
    mu(gid) = vals(id);
    if hasVars && strcmpi(opts.binMean,'mean')
        stdv(gid) = sqrt(vars(id));
        sem(gid) = sqrt(vars(id));
    else
        stdv(gid) = 0;
        sem(gid) = NaN;
    end
end

multiIdx = find(validOcc & nOcc > 1);
for jj = 1:numel(multiIdx)
    j = multiIdx(jj);
    g = groupIDs(j);
    id = ord(starts(j):stops(j));
    if hasVars
        vg = vars(id);
    else
        vg = [];
    end
    [mu(g), sem(g), stdv(g)] = local_bin_stats(vals(id), vg, opts);
end
end

% ======================================================================
function sem = local_apply_sem_floor(sem, np, opts)
good = isfinite(sem) & np >= opts.minPixPerBin;
if ~any(good)
    return
end
semRef = median(sem(good), 'omitnan');
if isfinite(semRef)
    semFloor = max(opts.semFloorAbs, opts.semFloorFrac*semRef);
else
    semFloor = opts.semFloorAbs;
end
if isfinite(semFloor) && semFloor > 0
    sem(good) = max(sem(good), semFloor);
end
end

% ======================================================================
function mode2 = local_empty_mode2(qCenters, qEdges)
nQ = numel(qCenters);
mode2 = struct();
mode2.qEdges = qEdges;
mode2.qCenters = qCenters(:);
mode2.phiEdges = cell(nQ,1);
mode2.phiCenters = cell(nQ,1);
mode2.Imean = cell(nQ,1);
mode2.Isem = cell(nQ,1);
mode2.Istd = cell(nQ,1);
mode2.Npix = cell(nQ,1);
mode2.wPhi = cell(nQ,1);
mode2.Nphi = nan(nQ,1);
mode2.Iq = nan(nQ,1);
mode2.IqSem = nan(nQ,1);
mode2.IqStd = nan(nQ,1);
mode2.Nq = zeros(nQ,1);
end

% ======================================================================
function mode1 = local_empty_mode1(qCenters, qEdges, phiEdges)
phiCenters = phiEdges(1:end-1) + 0.5*diff(phiEdges);
nQ = numel(qCenters);
nPhi = numel(phiCenters);
mode1 = struct();
mode1.qEdges = qEdges;
mode1.qCenters = qCenters(:);
mode1.phiEdges = phiEdges(:);
mode1.phiCenters = phiCenters(:);
mode1.Iqp = nan(nQ,nPhi);
mode1.Semqp = nan(nQ,nPhi);
mode1.Stdqp = nan(nQ,nPhi);
mode1.Nqp = zeros(nQ,nPhi);
mode1.Wqp = nan(nQ,nPhi);
end

% ======================================================================
function fitReady = local_make_fit_ready(mode2, opts)
nQ = numel(mode2.qCenters);
fitReady = struct();
fitReady.qCenters = mode2.qCenters;
fitReady.phi = cell(nQ,1);
fitReady.y = cell(nQ,1);
fitReady.sig = cell(nQ,1);
fitReady.w = cell(nQ,1);
fitReady.Npix = cell(nQ,1);
for iq = 1:nQ
    phi = mode2.phiCenters{iq};
    y = mode2.Imean{iq};
    s = mode2.Isem{iq};
    w = mode2.wPhi{iq};
    np = mode2.Npix{iq};
    if isempty(phi)
        continue
    end
    good = isfinite(phi) & isfinite(y) & isfinite(s) & isfinite(w) & ...
        w > 0 & np >= opts.minPixPerBin;
    fitReady.phi{iq} = phi(good);
    fitReady.y{iq} = y(good);
    fitReady.sig{iq} = s(good);
    fitReady.w{iq} = w(good);
    fitReady.Npix{iq} = np(good);
end
end

% ======================================================================
function [mu, sem, stdv] = local_bin_stats(vals, vars, opts)
vals = vals(:);
hasVars = ~isempty(vars);
if hasVars
    vars = vars(:);
    if numel(vars) ~= numel(vals)
        error('Variance and intensity vectors must have the same length.');
    end
end
goodVals = isfinite(vals);
vals = vals(goodVals);
if hasVars
    vars = vars(goodVals);
end
n = numel(vals);
mu = NaN;
sem = NaN;
stdv = NaN;
if n == 0
    return
end
useExactVarPropagation = hasVars && strcmpi(opts.binMean,'mean');

switch lower(opts.binMean)
    case 'mean'
        mu = mean(vals,'omitnan');
    case 'trimmean'
        if n >= 5
            mu = local_trimmean(vals, opts.trimPct);
        else
            mu = mean(vals,'omitnan');
        end
    otherwise
        error('opts.binMean must be mean or trimmean.');
end

if useExactVarPropagation
    sem = sqrt(sum(vars,'omitnan')) / n;
    stdv = sqrt(mean(vars,'omitnan'));
    return
end
if n == 1
    stdv = 0;
    sem = NaN;
    return
end

switch lower(opts.errModel)
    case 'stdsem'
        stdv = std(vals,0,'omitnan');
    case 'madsem'
        medv = median(vals,'omitnan');
        stdv = 1.4826*median(abs(vals-medv),'omitnan');
    otherwise
        error('opts.errModel must be stdsem or madsem.');
end
sem = stdv/sqrt(n);
end

% ======================================================================
function Nphi = local_choose_nphi(nPix, opts)
Nphi = round(nPix/max(opts.targetPixPerBin,1));
Nphi = max(Nphi,opts.minPhiBins);
Nphi = min(Nphi,opts.maxPhiBins);
NphiMax = max(1,floor(nPix/max(opts.minPixPerBin,1)));
Nphi = min(Nphi,NphiMax);
if opts.forceEvenPhiBins && Nphi > 1 && mod(Nphi,2) ~= 0
    Nphi = min(Nphi+1,NphiMax);
    if mod(Nphi,2) ~= 0 && Nphi > 2
        Nphi = Nphi-1;
    end
end
Nphi = max(Nphi,1);
end

% ======================================================================
function ang = local_wrap_to_center(ang, phiCtr)
ang = phiCtr + mod(ang-phiCtr+pi,2*pi)-pi;
end

% ======================================================================
function y = local_trimmean(x, pct)
%toolbox-independent symmetric trimmed mean.

x = x(isfinite(x));
if isempty(x)
    y = NaN;
    return
end

x = sort(x(:));
n = numel(x);
k = floor((pct/100) * n / 2);
if 2*k >= n
    y = mean(x);
else
    y = mean(x(k+1:n-k));
end
end

% ======================================================================
function out = local_merge_struct(def, in)
out = def;
if nargin < 2 || isempty(in)
    return
end
if ~isstruct(in)
    error('Configuration input must be a struct.');
end
names = fieldnames(in);
for k = 1:numel(names)
    name = names{k};
    if isstruct(in.(name)) && isfield(out,name) && isstruct(out.(name))
        out.(name) = local_merge_struct(out.(name),in.(name));
    else
        out.(name) = in.(name);
    end
end
end
