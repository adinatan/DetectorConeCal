function out = fit_pilatus100k_lut_robust(matFile, userCfg)
% Fast, accuracy-preserving PILATUS 100K calibrant (LaB6) pose fit.
% The expensive LUT stage is evaluated in vectorized candidate blocks on
% the CPU or GPU.  The coarse LUT may use a spatially stratified subset of
% ring pixels, but all promising candidates are rescored with all selected
% pixels in double precision before local refinement.  Final assignment and
% optimization always use the exact sorted LaB6 peak list and double
% precision.
%
% Basic use:
%   out = fit_pilatus100k_lut_robust('scan2_0000.mat');
%
% CPU/GPU control:
%   cfg.compute.device = 'auto'; % benchmark CPU/GPU and choose
%   cfg.compute.device = 'cpu';
%   cfg.compute.device = 'gpu';
%
% The output parameter convention is unchanged:
%   p = [beamX_px beamY_px zBeam_mm thetaX_deg thetaY_deg]
%
% Required dependency for the CIF reference:
%   cif2powder_CM.m
%
% Optional dependencies:
%   tripol_qphi_fast.m          fast q-phi output
%   an_detector_core.m         fallback q-phi output
%
% Adi Natan , 2026 SSRL beamtime.

if nargin < 1 || isempty(matFile)
    matFile = 'scan2_0000.mat';
end
if nargin < 2 || isempty(userCfg)
    userCfg = struct();
end

tTotal = tic;
timing = struct();

cfg = local_default_cfg();
cfg = local_merge_struct(cfg, userCfg);
rng(cfg.rngSeed);

thisDir = fileparts(mfilename('fullpath'));
if ~isempty(thisDir)
    addpath(thisDir);
end

computeInfo = local_resolve_compute(cfg.compute);

tPhase = tic;
S = load(matFile);
[Im0, imName] = local_find_image(S);
Im0 = double(Im0);
Im0(~isfinite(Im0)) = 0;
[Im, orientInfo] = local_apply_orientation(Im0, cfg.image.orientation);
det = local_detector_from_image(size(Im), cfg);
timing.load = toc(tPhase);

tPhase = tic;
ref = local_make_reference(cfg);
timing.reference = toc(tPhase);

tPhase = tic;
prep = local_prepare_image(Im, cfg);
ptsFull = local_make_score_points(prep, cfg);
ptsCoarse = local_make_coarse_score_points(ptsFull, cfg.lut.coarsePointCount);
timing.preprocessing = toc(tPhase);

fprintf('\n=== PILATUS 100K fast old+LUT pose fit ===\n');
fprintf('MAT file: %s, variable: %s\n', matFile, imName);
fprintf('Image orientation: %s; image used for fit: %d rows x %d cols\n', ...
    orientInfo.name, size(Im,1), size(Im,2));
fprintf('Pixel size: %.6g mm; z range %.1f--%.1f mm; prior tilt %.2f deg about %s\n', ...
    det.pixel_mm, cfg.bounds.zBeam_mm(1), cfg.bounds.zBeam_mm(2), ...
    cfg.prior.tiltDeg, cfg.prior.mainTiltAxis);
fprintf('Signal pixels: %d full, %d coarse; threshold %.4g (sigma %.2f, percentile %.2f, relaxed=%d)\n', ...
    numel(ptsFull.x), numel(ptsCoarse.x), prep.threshold, ...
    prep.thresholdSigmaSetting, prep.thresholdPercentileSetting, prep.thresholdWasRelaxed);
fprintf('Reference: %s; %d LaB6 Q peaks %.3f--%.3f A^-1\n', ...
    ref.sourceText, numel(ref.qPeaks_Ainv), min(ref.qPeaks_Ainv), max(ref.qPeaks_Ainv));
fprintf('Requested compute device: %s; GPU available: %d\n', ...
    computeInfo.requestedDevice, computeInfo.gpuAvailable);

cfgRequested = cfg;
[coarse, fit, computeInfo, searchHistory, cfg] = local_robust_search( ...
    ptsCoarse, ptsFull, ref, det, cfg, computeInfo);
timing.lut = sum([searchHistory.lutTime_s]);
timing.refinement = sum([searchHistory.refineTime_s]);

pFit = fit.bestP;

tPhase = tic;
geom = local_build_geom_old(size(Im), det.pixel_mm, pFit, cfg.beam.EkeV, cfg.beam.polMode);
profileRaw = local_profile_from_geom(Im, geom, cfg.profile.qEdges_Ainv, cfg.profile.trimPct);
profileEnh = local_profile_from_geom(prep.enhancedForProfile, geom, cfg.profile.qEdges_Ainv, cfg.profile.trimPct);
pose6d = local_old_params_to_pose6d(pFit, size(Im), det.pixel_mm);
timing.outputMaps = toc(tPhase);

tripCfg = local_make_tripol_cfg(cfg, prep.validMask);
TP = [];
timing.tripol = 0;
if cfg.output.runTripol
    tPhase = tic;
    engine = lower(char(cfg.output.tripolEngine));
    try
        switch engine
            case 'fast'
                if exist('tripol_qphi_fast', 'file') == 2
                    TP = tripol_qphi_fast(Im, geom, tripCfg);
                    fprintf('Built q-phi output with tripol_qphi_fast.\n');
                elseif exist('an_detector_core', 'file') == 2
                    warning('tripol_qphi_fast not found; using an_detector_core instead.');
                    TP = an_detector_core('tripol', Im, geom, tripCfg);
                else
                    warning('No q-phi implementation was found; out.TP is empty.');
                end
            case 'aps'
                if exist('an_detector_core', 'file') ~= 2
                    error('an_detector_core.m was not found.');
                end
                TP = an_detector_core('tripol', Im, geom, tripCfg);
                fprintf('Built q-phi output with an_detector_core.\n');
            case {'none','off'}
                % Explicitly disabled.
            otherwise
                error('cfg.output.tripolEngine must be fast, aps, or none.');
        end
    catch ME
        warning('Could not build q-phi output: %s', ME.message);
        TP = [];
    end
    timing.tripol = toc(tPhase);
end

out = struct();
out.matFile = matFile;
out.imageVariable = imName;
out.inputImage = Im0;
out.imageUsed = Im;
out.orientation = orientInfo;
out.cfg = cfg;
out.cfgRequested = cfgRequested;
out.searchHistory = searchHistory;
out.compute = computeInfo;
out.detector = det;
out.ref = ref;
out.prep = prep;
out.scorePoints = ptsFull;
out.coarseScorePoints = ptsCoarse;
out.coarse = coarse;
out.fit = fit;
out.params = pFit;
out.params_columns = {'beamX_px','beamY_px','zBeam_mm','theta_x_deg','theta_y_deg'};
out.pose6d = pose6d;
out.pose6d_columns = {'Tx_mm','Ty_mm','Tz_mm','theta_x_deg','theta_y_deg','theta_z_deg'};
out.geom = geom;
out.profileRaw = profileRaw;
out.profileEnhanced = profileEnh;
out.tripCfg = tripCfg;
out.TP = TP;
out.validation = local_validate_result(pFit, fit, cfg.validation);

fprintf('\nBest fitted old-geometry parameters:\n');
fprintf('  p = [beamX beamY z tx ty] = [%+.4f %+.4f %+.4f  %+.4f %+.4f]\n', pFit);
fprintf('  pose6d detector center [Tx Ty Tz rx ry rz] = [%+.4f %+.4f %+.4f  %+.4f %+.4f %+.4f]\n', pose6d);
fprintf('  final score = %.6f; raw peak score = %.6f; q residual WRMS = %.5f A^-1; assigned pixels = %d\n', ...
    fit.bestScore, fit.bestPeakScore, fit.bestAssignment.wrms_Ainv, nnz(fit.bestAssignment.keep));
fprintf('  LUT device = %s; search rounds = %d\n', computeInfo.activeDevice, numel(searchHistory));

if cfg.output.makePlots
    tPhase = tic;
    local_plot_signal_processing(out);
    local_plot_pose_result(out);
    timing.plots = toc(tPhase);
else
    timing.plots = 0;
end

timing.totalBeforeSave = toc(tTotal);
out.timing = timing;
out.timing.save = NaN;
out.timing.total = timing.totalBeforeSave;

if cfg.output.saveMat
    tPhase = tic;
    local_save_output(out, cfg.output);
    saveTime = toc(tPhase);
    fprintf('Saved %s in %.3f s\n', cfg.output.matFile, saveTime);
else
    saveTime = 0;
end

out.timing.save = saveTime;
out.timing.total = toc(tTotal);

fprintf('\nTiming summary (s): reference %.3f, prep %.3f, LUT %.3f, refine %.3f, maps %.3f, tripol %.3f, plots %.3f, save %.3f, total %.3f\n', ...
    out.timing.reference, out.timing.preprocessing, out.timing.lut, out.timing.refinement, ...
    out.timing.outputMaps, out.timing.tripol, out.timing.plots, out.timing.save, out.timing.total);
end
% ======================================================================
function cfg = local_default_cfg()
cfg = struct();
cfg.rngSeed = 11;

cfg.image = struct();
cfg.image.orientation = 'as_loaded';

cfg.detector = struct();
cfg.detector.name = 'PILATUS_100K';
cfg.detector.pixel_um = 172;

cfg.beam = struct();
cfg.beam.EkeV = 16.6;
cfg.beam.polMode = 'horizontal';

cfg.ref = struct();
cfg.ref.cifFile = '1000057(1).cif';
cfg.ref.altCifFile = '1000057.cif';
cfg.ref.qGrid_Ainv = 0.3:0.0025:9.5;
cfg.ref.minRelPeakI = 0.004;
cfg.ref.profileFwhmQ_Ainv = 0.035;
cfg.ref.useCache = true;
cfg.ref.cacheFile = 'lab6_reference_16p6keV_cache.mat';

cfg.prep = struct();
cfg.prep.capHighPercentile = 99.98;
cfg.prep.backgroundWindowPx = 35;
cfg.prep.thresholdSigma = 2.5;
cfg.prep.thresholdPercentile = 98.5;
cfg.prep.minObjectPixels = 2;
cfg.prep.useClean = true;
cfg.prep.maxScorePixels = 6000;
cfg.prep.weightOffset = 0.25;
cfg.prep.computeLegacyAudit = false;
% Robust mode can gently lower the two thresholds when a high-angle image
% leaves too few ring pixels.  The requested thresholds remain the starting
% point; they are never made more aggressive automatically.
cfg.prep.autoRelaxThreshold = true;
cfg.prep.minScorePixels = 1100;
cfg.prep.relaxSigmaFloor = 2.2;
cfg.prep.relaxPercentileFloor = 98.2;
cfg.prep.relaxSigmaStep = 0.15;
cfg.prep.relaxPercentileStep = 0.20;

cfg.bounds = struct();
cfg.bounds.beamX_px = [-60 95];
cfg.bounds.beamY_px = [35 155];
cfg.bounds.zBeam_mm = [140 190];
cfg.bounds.thetaX_deg = [-8 8];
cfg.bounds.thetaY_deg = [6 22];

cfg.prior = struct();
cfg.prior.mainTiltAxis = 'y';
cfg.prior.tiltDeg = 14;
cfg.prior.tiltSign = +1;
% In robust mode the user bounds are the hard bounds.  The prior guides the
% sampling and contributes only a bounded soft penalty unless this flag is
% explicitly enabled.
cfg.prior.constrainBounds = false;
cfg.prior.tiltSigma_deg = 12;
cfg.prior.tiltPenaltyWeight = 0.004;
cfg.prior.penaltyType = 'soft';       % 'soft', 'quadratic', or 'none'
cfg.prior.maxPenalty = 0.035;         % prevents a wrong prior hiding a good basin

cfg.compute = struct();
cfg.compute.device = 'auto';              % 'auto', 'cpu', or 'gpu'
cfg.compute.autoBenchmark = true;
cfg.compute.autoBenchmarkCandidates = 2048;
cfg.compute.autoGpuMinSpeedup = 1.05;
cfg.compute.cpuBlockCandidates = 512;
cfg.compute.gpuBlockCandidates = 4096;
cfg.compute.gpuPrecision = 'single';       % LUT only; final fit is double
cfg.compute.cpuParallel = 'auto';          % false, true, or 'auto' (existing pool only)
cfg.compute.hitLookupDQ_Ainv = 1e-5;
cfg.compute.gpuRequired = false;
cfg.compute.resetGPU = false;

cfg.lut = struct();
cfg.lut.nGlobal = 90000;
% Multiple local levels are substantially more reliable than one very narrow
% cloud after a broad 5-D LUT.  Scalars are accepted and are repeated.
cfg.lut.nLocalSeeds = [8 5 3];
cfg.lut.nLocalPerSeed = [3000 2200 1400];
cfg.lut.topK = 40;
cfg.lut.localHalfWidth = [ ...
    45 32 28 12 14; ...
    18 14 12  5  6; ...
     6  6  4  2  2];
cfg.lut.peakSigmaQ_Ainv = 0.035;
cfg.lut.progressEvery = 10000;
cfg.lut.coarsePointCount = 450;
cfg.lut.globalFullRescoreTopK = 1800;
cfg.lut.localFullRescoreTopK = [1200 800 500];
cfg.lut.seedMinDistance = [0.80 0.65 0.50];
cfg.lut.useGridIfSmall = true;
cfg.lut.gridStep = [8 8 2 4 2];
cfg.lut.maxGridCandidates = 60000;
% Compact deterministic prior grid around the known y-tilt basin.
cfg.lut.priorGridXCount = 9;
cfg.lut.priorGridYCount = 5;
cfg.lut.priorGridZStep_mm = 10;
cfg.lut.priorMainTiltOffsets_deg = -10:5:10;
cfg.lut.priorMinorTiltOffsets_deg = [-4 0 4];
cfg.lut.includeBroadAnchorGrid = true;
cfg.lut.broadAnchorXCount = 5;
cfg.lut.broadAnchorYCount = 3;
cfg.lut.broadAnchorZCount = 5;
cfg.lut.broadAnchorTiltStep_deg = 15;
cfg.lut.sobolSkipOffset = 0;

cfg.assign = struct();
cfg.assign.maxDQ_Ainv = 0.075;
cfg.assign.minAssignedPixels = 80;
cfg.assign.maxRefinePixels = 5000;

cfg.opt = struct();
cfg.opt.topKToRefine = 3;
cfg.opt.refineMinDistance = 0.18;
cfg.opt.reassignCycles = 2;
cfg.opt.reassignParamTol = 2e-4;
cfg.opt.robustScaleQ_Ainv = 0.010;
cfg.opt.maxIter = 160;
cfg.opt.maxFunEvals = 900;
cfg.opt.solver = 'auto';                   % auto, lsqnonlin, fmincon, fminsearch
cfg.opt.stepTolerance = 1e-6;
cfg.opt.functionTolerance = 1e-8;
cfg.opt.optimalityTolerance = 1e-6;

cfg.search = struct();
cfg.search.autoExpandBounds = true;
cfg.search.maxRounds = 3;
cfg.search.enforceMinimumGlobalCandidates = true;
cfg.search.minimumGlobalCandidates = 60000;
cfg.search.edgeFraction = 0.035;
cfg.search.hardEdgeFraction = 0.010;
cfg.search.minAcceptablePeakScore = 0.82;
cfg.search.maxAcceptableWRMS_Ainv = 0.015;
cfg.search.resampleOnLowScore = true;
cfg.search.samplingGrowth = 1.8;
% Expansion is a fraction of the current span, separately for each side.
cfg.search.expandFraction = [1.0 0.60 0.60 0.50 0.50];
cfg.search.expandParameters = [true true true true true];
% Safety limits, rows are [lower upper] for [beamX beamY z tx ty].
cfg.search.absoluteLimits = [ ...
   -2500 2500; ...
   -1200 1200; ...
      20 1000; ...
     -82   82; ...
     -82   82];

cfg.profile = struct();
cfg.profile.qEdges_Ainv = 0.3:0.015:8.5;
cfg.profile.trimPct = 37;

cfg.tripol = struct();
cfg.tripol.qEdges = 0.3:0.025:8.5;
cfg.tripol.fixedPhiEdges = linspace(-pi, pi, 361);
cfg.tripol.applyTotalCorr = true;
cfg.tripol.minPixPerBin = 1;
cfg.tripol.minPixPerRing = 1;
cfg.tripol.targetPixPerBin = 80;
cfg.tripol.binMean = 'trimmean';
cfg.tripol.trimPct = 37;

cfg.validation = struct();
cfg.validation.referenceP = [];
cfg.validation.tolerance = [0.10 0.10 0.05 0.05 0.05];
cfg.validation.warnOnly = true;

cfg.output = struct();
cfg.output.makePlots = true;
cfg.output.saveMat = true;
cfg.output.matFile = 'scan2_pilatus100k_oldlut_fast_pose_fit.mat';
cfg.output.saveVersion = '-v7.3';
cfg.output.saveCompact = false;
cfg.output.runTripol = true;
cfg.output.tripolEngine = 'fast';
cfg.output.nContourPeaks = 12;
end

% ======================================================================
function det = local_detector_from_image(imSize, cfg)
det = struct();
det.name = cfg.detector.name;
det.nRows = imSize(1);
det.nCols = imSize(2);
det.pixel_um = cfg.detector.pixel_um;
det.pixel_mm = det.pixel_um / 1000;
det.activeSize_mm = [det.nCols det.nRows] * det.pixel_mm;
det.center_px = [(det.nCols+1)/2, (det.nRows+1)/2];
end

% ======================================================================
function [Im, info] = local_apply_orientation(Im0, orientation)
name = lower(char(orientation));
switch name
    case {'as_loaded','none','raw'}
        Im = Im0;
        outName = 'as_loaded';
        transform = 'identity';
    case {'long_axis_vertical','vertical','lab_y_long'}
        if size(Im0,2) > size(Im0,1)
            Im = rot90(Im0, 1);
            transform = 'rot90ccw because columns were the long axis';
        else
            Im = Im0;
            transform = 'identity because rows were already the long axis';
        end
        outName = 'long_axis_vertical';
    case {'rot90ccw','ccw'}
        Im = rot90(Im0, 1);
        outName = 'rot90ccw';
        transform = 'rot90ccw';
    case {'rot90cw','cw'}
        Im = rot90(Im0, -1);
        outName = 'rot90cw';
        transform = 'rot90cw';
    case {'transpose'}
        Im = Im0.';
        outName = 'transpose';
        transform = 'transpose';
    otherwise
        error('Unknown cfg.image.orientation: %s', orientation);
end
info = struct('name', outName, 'transform', transform, 'inputSize', size(Im0), 'fitSize', size(Im));
end

% ======================================================================
function [Im, name] = local_find_image(S)
if isfield(S, 'Im') && isnumeric(S.Im) && ismatrix(S.Im)
    Im = S.Im;
    name = 'Im';
    return
end
names = fieldnames(S);
for i = 1:numel(names)
    v = S.(names{i});
    if isnumeric(v) && ismatrix(v) && numel(v) > 1000
        Im = v;
        name = names{i};
        return
    end
end
error('Could not find a numeric 2D detector image in the MAT file.');
end

% ======================================================================
function prep = local_prepare_image(Im, cfg)
I0 = double(Im);
valid = isfinite(I0) & I0 > 0;
if nnz(valid) < 1000
    error('Too few positive/valid detector pixels.');
end

cap = local_prctile(I0(valid), cfg.prep.capHighPercentile);
Iclip = min(max(I0, 0), cap);
L = log1p(Iclip);
fillValue = median(L(valid));
Lfill = L;
Lfill(~valid | ~isfinite(Lfill)) = fillValue;

w = max(5, round(cfg.prep.backgroundWindowPx));
if mod(w,2) == 0
    w = w + 1;
end
bg = local_movmedian2(Lfill, w);
E = Lfill - bg;
E(~valid) = NaN;
Epos = max(E, 0);
vals = Epos(valid & isfinite(Epos));
medv = median(vals);
sigv = 1.4826 * median(abs(vals - medv));
if ~isfinite(sigv) || sigv <= 0
    sigv = std(vals);
end
sigmaSetting = cfg.prep.thresholdSigma;
pctSetting = cfg.prep.thresholdPercentile;
[thrSigma, thrPct, thr, BWraw, BW] = local_threshold_ring_mask( ...
    Epos, valid, vals, medv, sigv, sigmaSetting, pctSetting, cfg.prep);

if cfg.prep.autoRelaxThreshold
    targetN = min(cfg.prep.minScorePixels, cfg.prep.maxScorePixels);
    while nnz(BW) < targetN && ...
            (sigmaSetting > cfg.prep.relaxSigmaFloor || ...
             pctSetting > cfg.prep.relaxPercentileFloor)
        sigmaSetting = max(cfg.prep.relaxSigmaFloor, ...
            sigmaSetting - cfg.prep.relaxSigmaStep);
        pctSetting = max(cfg.prep.relaxPercentileFloor, ...
            pctSetting - cfg.prep.relaxPercentileStep);
        [thrSigma, thrPct, thr, BWraw, BW] = local_threshold_ring_mask( ...
            Epos, valid, vals, medv, sigv, sigmaSetting, pctSetting, cfg.prep);
    end
end

% This legacy view is diagnostic only and is skipped by default.
im3old = [];
if cfg.prep.computeLegacyAudit
    im3old = zeros(size(I0));
    try
        CC = bwlabel(valid);
        for k = 1:max(CC(:))
            id = CC == k;
            base = local_trimmean(I0(id), 37);
            im3old(id) = I0(id) - base;
        end
    catch
        base = local_trimmean(I0(valid), 37);
        im3old(valid) = I0(valid) - base;
    end
    im3old = max(im3old, 0);
end

prep = struct();
prep.validMask = valid;
prep.raw = I0;
prep.clipped = Iclip;
prep.logImage = L;
prep.background = bg;
prep.enhanced = Epos;
prep.enhancedForProfile = Epos;
prep.oldHighpass = im3old;
prep.ringMaskRaw = BWraw;
prep.ringMask = BW;
prep.threshold = thr;
prep.thresholdSigma = thrSigma;
prep.thresholdPercentile = thrPct;
prep.thresholdSigmaSetting = sigmaSetting;
prep.thresholdPercentileSetting = pctSetting;
prep.thresholdWasRelaxed = sigmaSetting < cfg.prep.thresholdSigma || ...
    pctSetting < cfg.prep.thresholdPercentile;
prep.backgroundWindowPx = w;
end

% ======================================================================
function [thrSigma, thrPct, thr, BWraw, BW] = local_threshold_ring_mask( ...
    Epos, valid, vals, medv, sigv, sigmaSetting, pctSetting, prepCfg)
thrSigma = medv + sigmaSetting * sigv;
thrPct = local_prctile(vals, pctSetting);
thr = max(thrSigma, thrPct);
BWraw = valid & isfinite(Epos) & Epos >= thr;
BW = BWraw;
if prepCfg.useClean
    if exist('bwareaopen','file') == 2
        BW = bwareaopen(BW, prepCfg.minObjectPixels);
    end
    if exist('bwmorph','file') == 2
        BW = bwmorph(BW, 'clean');
    end
end
end

% ======================================================================
function bg = local_movmedian2(A, w)
try
    bg = movmedian(A, w, 2, 'omitnan');
    bg = movmedian(bg, w, 1, 'omitnan');
catch
    bg = movmedian(A, w, 2);
    bg = movmedian(bg, w, 1);
end
end

% ======================================================================
function pts = local_make_score_points(prep, cfg)
[y, x] = find(prep.ringMask);
if numel(x) < 50
    error('Only %d ring-enhanced pixels found. Lower cfg.prep.thresholdPercentile or cfg.prep.thresholdSigma.', numel(x));
end
x = double(x(:));
y = double(y(:));
w = prep.enhanced(prep.ringMask);
w = double(w(:));
w(~isfinite(w)) = 0;
if max(w) > 0
    w = w / max(w);
else
    w = ones(size(w));
end
w = cfg.prep.weightOffset + (1-cfg.prep.weightOffset) * w;

if numel(x) > cfg.prep.maxScorePixels
    [~, ord] = sortrows([x y]);
    take = unique(round(linspace(1, numel(ord), cfg.prep.maxScorePixels)));
    idx = ord(take);
    x = x(idx);
    y = y(idx);
    w = w(idx);
end
pts = struct('x', x, 'y', y, 'w', w);
end

% ======================================================================
function ptsOut = local_make_coarse_score_points(ptsIn, nWanted)
n = numel(ptsIn.x);
if isempty(nWanted) || nWanted <= 0 || n <= nWanted
    ptsOut = ptsIn;
    return
end

% Deterministic spatially stratified selection.  Sorting by x and then y
% retains samples from every visible arc rather than simply keeping the
% highest-intensity pixels.
[~, ord] = sortrows([ptsIn.x(:), ptsIn.y(:)], [1 2]);
take = unique(round(linspace(1, n, nWanted)));
idx = ord(take);
ptsOut = struct('x', ptsIn.x(idx), 'y', ptsIn.y(idx), 'w', ptsIn.w(idx));
end

% ======================================================================
function ref = local_make_reference(cfg)
qGrid = cfg.ref.qGrid_Ainv(:);
EkeV = cfg.beam.EkeV;
lambda_A = 12.398419843320026 / EkeV;
qPeaks = [];
peakI = [];
sourceText = 'fallback cubic LaB6 Q peak list';

cifFile = local_resolve_file(cfg.ref.cifFile, cfg.ref.altCifFile);
cacheKey = local_reference_cache_key(cifFile, cfg);
if cfg.ref.useCache && ~isempty(cfg.ref.cacheFile) && exist(cfg.ref.cacheFile, 'file') == 2
    try
        C = load(cfg.ref.cacheFile, 'refCache');
        if isfield(C, 'refCache') && isfield(C.refCache, 'key') && ...
                isequaln(C.refCache.key, cacheKey)
            ref = C.refCache.ref;
            ref.sourceText = [ref.sourceText ' [cached]'];
            return
        end
    catch ME
        warning('Could not read LaB6 reference cache: %s', ME.message);
    end
end

if ~isempty(cifFile) && exist('cif2powder_CM', 'file') == 2
    try
        opts = struct();
        opts.twoThetaGrid_deg = 0.05:0.02:120;
        opts.qGrid_Ainv = qGrid(:).';
        opts.profile.widthModel = 'fixed';
        opts.profile.fixedFWHM_deg = 0.08;
        opts.profile.shape = 'pvoigt';
        opts.profile.eta = 0.35;
        opts.corrections.useLorentz = true;
        opts.corrections.usePolarization = false;
        opts.verbose = false;
        powder = cif2powder_CM(cifFile, EkeV, opts);
        qPeaks = powder.linesMerged.Q_Ainv(:);
        peakI = powder.linesMerged.I(:);
        sourceText = sprintf('cif2powder_CM(%s, %.3g keV)', cifFile, EkeV);
    catch ME
        warning('Could not build LaB6 reference from CIF: %s. Using fallback cubic list.', ME.message);
    end
end

if isempty(qPeaks)
    a_A = 4.1570;
    hmax = 14;
    nList = [];
    for h = 0:hmax
        for k = 0:hmax
            for l = 0:hmax
                nval = h*h + k*k + l*l;
                if nval > 0
                    nList(end+1,1) = nval;  
                end
            end
        end
    end
    nList = unique(nList);
    qPeaks = (2*pi/a_A) * sqrt(nList(:));
    peakI = 1 ./ (1 + (qPeaks/4).^2);
end

arg = qPeaks * lambda_A / (4*pi);
keep = isfinite(qPeaks) & isfinite(peakI) & peakI > 0 & arg > 0 & arg < 1 & ...
    qPeaks >= min(qGrid) & qPeaks <= max(qGrid);
qPeaks = qPeaks(keep);
peakI = peakI(keep);
peakI = peakI / max(peakI + eps);
keep = peakI >= cfg.ref.minRelPeakI;
qPeaks = qPeaks(keep);
peakI = peakI(keep);
[qPeaks, ord] = sort(qPeaks(:));
peakI = peakI(ord);

profile = zeros(size(qGrid));
for k = 1:numel(qPeaks)
    profile = profile + peakI(k) * local_gaussian_fwhm(qGrid, qPeaks(k), cfg.ref.profileFwhmQ_Ainv);
end
profile = profile / max(profile + eps);

ref = struct();
ref.lambda_A = lambda_A;
ref.qGrid_Ainv = qGrid;
ref.qPeaks_Ainv = qPeaks;
ref.peakI = peakI;
ref.profile = profile(:);
ref.sourceText = sourceText;
ref.qMidEdges_Ainv = [-Inf; 0.5*(qPeaks(1:end-1) + qPeaks(2:end)); Inf];

if cfg.ref.useCache && ~isempty(cfg.ref.cacheFile)
    try
        refCache = struct('key', cacheKey, 'ref', ref);   
        save(cfg.ref.cacheFile, 'refCache', '-v7');
    catch ME
        warning('Could not save LaB6 reference cache: %s', ME.message);
    end
end
end

% ======================================================================
function key = local_reference_cache_key(cifFile, cfg)
key = struct();
key.cifFile = cifFile;
key.cifBytes = NaN;
key.cifDatenum = NaN;
if ~isempty(cifFile) && exist(cifFile, 'file') == 2
    d = dir(cifFile);
    key.cifBytes = d.bytes;
    key.cifDatenum = d.datenum;
end
key.EkeV = cfg.beam.EkeV;
key.qGrid_Ainv = cfg.ref.qGrid_Ainv(:).';
key.minRelPeakI = cfg.ref.minRelPeakI;
key.profileFwhmQ_Ainv = cfg.ref.profileFwhmQ_Ainv;
end

% ======================================================================
function [coarseBest, fitBest, computeInfo, history, cfgBest] = local_robust_search( ...
    ptsCoarse, ptsFull, ref, det, cfg, computeInfo)

cfgRound = cfg;
if cfg.search.enforceMinimumGlobalCandidates
    cfgRound.lut.nGlobal = max(cfgRound.lut.nGlobal, ...
        round(cfg.search.minimumGlobalCandidates));
end
history = repmat(struct( ...
    'round',[], 'bounds',[], 'nGlobal',[], 'lutTime_s',[], ...
    'refineTime_s',[], 'bestP',[], 'bestScore',[], 'peakScore',[], ...
    'wrms_Ainv',[], 'nAssigned',[], 'action','', 'reason',''), 0, 1);

coarseBest = [];
fitBest = [];
cfgBest = cfg;
bestMetric = -Inf;

nRounds = 1;
if cfg.search.autoExpandBounds || cfg.search.resampleOnLowScore
    nRounds = max(1, round(cfg.search.maxRounds));
end

for iRound = 1:nRounds
    [lbRound, ubRound] = local_global_bounds(cfgRound);
    fprintf('\n============================================================\n');
    fprintf('Robust search round %d/%d\n', iRound, nRounds);
    fprintf('  bounds lower = [%s]\n', num2str(lbRound,' %.5g'));
    fprintf('  bounds upper = [%s]\n', num2str(ubRound,' %.5g'));
    fprintf('  nGlobal = %d; prior hard-bound constraint = %d\n', ...
        cfgRound.lut.nGlobal, logical(cfgRound.prior.constrainBounds));

    t0 = tic;
    [coarseRound, computeInfo] = local_lut_search( ...
        ptsCoarse, ptsFull, ref, det, cfgRound, computeInfo);
    lutTime = toc(t0);

    t0 = tic;
    fitRound = local_refine_from_lut(coarseRound, ptsFull, ref, det, cfgRound);
    refineTime = toc(t0);

    wrms = fitRound.bestAssignment.wrms_Ainv;
    peakScore = fitRound.bestPeakScore;
    metric = fitRound.bestScore;

    h = struct();
    h.round = iRound;
    h.bounds = [coarseRound.lb(:), coarseRound.ub(:)];
    h.nGlobal = cfgRound.lut.nGlobal;
    h.lutTime_s = lutTime;
    h.refineTime_s = refineTime;
    h.bestP = fitRound.bestP;
    h.bestScore = fitRound.bestScore;
    h.peakScore = peakScore;
    h.wrms_Ainv = wrms;
    h.nAssigned = nnz(fitRound.bestAssignment.keep);
    h.action = '';
    h.reason = '';

    if metric > bestMetric || isempty(fitBest)
        bestMetric = metric;
        coarseBest = coarseRound;
        fitBest = fitRound;
        cfgBest = cfgRound;
    end

    [continueSearch, cfgNext, action, reason] = local_next_search_round( ...
        fitRound, coarseRound, cfgRound, iRound, nRounds);
    h.action = action;
    h.reason = reason;
    history(end+1,1) = h;  

    fprintf('Round %d result: peak score %.6f, WRMS %.5f A^-1, p=[%s]\n', ...
        iRound, peakScore, wrms, num2str(fitRound.bestP,' %+.4f'));
    if ~isempty(action)
        fprintf('Next action: %s. %s\n', action, reason);
    end

    if ~continueSearch
        break
    end
    cfgRound = cfgNext;
end

if isempty(fitBest)
    error('Robust search did not produce a valid fitted candidate.');
end
end

% ======================================================================
function [continueSearch, cfgNext, action, reason] = local_next_search_round( ...
    fit, coarse, cfg, iRound, nRounds)

continueSearch = false;
cfgNext = cfg;
action = '';
reason = '';

if iRound >= nRounds
    return
end

p = fit.bestP(:).';
lb = coarse.lb(:).';
ub = coarse.ub(:).';
span = max(ub-lb, eps);
fLo = (p-lb)./span;
fHi = (ub-p)./span;

peakLow = fit.bestPeakScore < cfg.search.minAcceptablePeakScore;
wrmsHigh = fit.bestAssignment.wrms_Ainv > cfg.search.maxAcceptableWRMS_Ainv;
qualityLow = peakLow || wrmsHigh;

nearLo = fLo <= cfg.search.edgeFraction;
nearHi = fHi <= cfg.search.edgeFraction;
hardLo = fLo <= cfg.search.hardEdgeFraction;
hardHi = fHi <= cfg.search.hardEdgeFraction;

expandLo = hardLo | (qualityLow & nearLo);
expandHi = hardHi | (qualityLow & nearHi);
canExpand = logical(cfg.search.expandParameters(:).');
expandLo = expandLo & canExpand;
expandHi = expandHi & canExpand;

if cfg.prior.constrainBounds
    % A hard prior is an explicit request.  Do not silently expand the angle
    % dimensions against it, although positional bounds may still expand.
    expandLo(4:5) = false;
    expandHi(4:5) = false;
end

if cfg.search.autoExpandBounds && any(expandLo | expandHi)
    frac = local_expand_to_five(cfg.search.expandFraction, 'cfg.search.expandFraction');
    absLim = double(cfg.search.absoluteLimits);
    if ~isequal(size(absLim),[5 2])
        error('cfg.search.absoluteLimits must be a 5-by-2 matrix.');
    end

    newLb = lb;
    newUb = ub;
    for j = 1:5
        if expandLo(j)
            newLb(j) = max(absLim(j,1), lb(j) - frac(j)*span(j));
        end
        if expandHi(j)
            newUb(j) = min(absLim(j,2), ub(j) + frac(j)*span(j));
        end
    end

    if any(newLb < lb) || any(newUb > ub)
        cfgNext = local_set_cfg_bounds(cfgNext, newLb, newUb);
        cfgNext.lut.sobolSkipOffset = cfg.lut.sobolSkipOffset + ...
            cfg.lut.nGlobal + 1000000*iRound;
        continueSearch = true;
        action = 'expand bounds';
        names = {'beamX','beamY','z','thetaX','thetaY'};
        parts = {};
        for j = 1:5
            if expandLo(j)
                parts{end+1} = sprintf('%s lower', names{j});  
            end
            if expandHi(j)
                parts{end+1} = sprintf('%s upper', names{j});  
            end
        end
        reason = sprintf('fit/seed evidence is near %s; qualityLow=%d', ...
            strjoin(parts, ', '), qualityLow);
        return
    end
end

if qualityLow && cfg.search.resampleOnLowScore
    cfgNext.lut.nGlobal = max(cfg.lut.nGlobal+1, ...
        ceil(cfg.lut.nGlobal * cfg.search.samplingGrowth));
    cfgNext.lut.globalFullRescoreTopK = max(cfg.lut.globalFullRescoreTopK, ...
        ceil(cfg.lut.globalFullRescoreTopK * sqrt(cfg.search.samplingGrowth)));
    cfgNext.lut.sobolSkipOffset = cfg.lut.sobolSkipOffset + ...
        cfg.lut.nGlobal + 1000000*iRound;
    % A poor fit should not be repeatedly pulled toward an incorrect guess.
    cfgNext.prior.tiltPenaltyWeight = 0.5 * cfg.prior.tiltPenaltyWeight;
    continueSearch = true;
    action = 'increase sampling';
    reason = sprintf('peakScore %.4f (target %.4f), WRMS %.5f (target %.5f)', ...
        fit.bestPeakScore, cfg.search.minAcceptablePeakScore, ...
        fit.bestAssignment.wrms_Ainv, cfg.search.maxAcceptableWRMS_Ainv);
end
end

% ======================================================================
function v = local_expand_to_five(vIn, name)
v = double(vIn(:).');
if isscalar(v)
    v = repmat(v,1,5);
elseif numel(v) ~= 5
    error('%s must be scalar or contain five values.', name);
end
end

% ======================================================================
function cfg = local_set_cfg_bounds(cfg, lb, ub)
cfg.bounds.beamX_px = [lb(1) ub(1)];
cfg.bounds.beamY_px = [lb(2) ub(2)];
cfg.bounds.zBeam_mm = [lb(3) ub(3)];
cfg.bounds.thetaX_deg = [lb(4) ub(4)];
cfg.bounds.thetaY_deg = [lb(5) ub(5)];
end

% ======================================================================
function [coarse, computeInfo] = local_lut_search(ptsCoarse, ptsFull, ref, det, cfg, computeInfo)
[lb, ub] = local_global_bounds(cfg);

Gprior = local_prior_grid_candidates(lb, ub, cfg);
Ggrid = zeros(0,5);
if cfg.lut.useGridIfSmall
    xVec = lb(1):cfg.lut.gridStep(1):ub(1);
    yVec = lb(2):cfg.lut.gridStep(2):ub(2);
    zVec = lb(3):cfg.lut.gridStep(3):ub(3);
    txVec = lb(4):cfg.lut.gridStep(4):ub(4);
    tyVec = lb(5):cfg.lut.gridStep(5):ub(5);
    nGrid = numel(xVec)*numel(yVec)*numel(zVec)*numel(txVec)*numel(tyVec);
    if nGrid <= cfg.lut.maxGridCandidates
        [X,Y,Z,TX,TY] = ndgrid(xVec,yVec,zVec,txVec,tyVec);
        Ggrid = [X(:) Y(:) Z(:) TX(:) TY(:)];
    else
        fprintf('Skipping full grid (%d candidates > maxGridCandidates %d); using anchor grid + Sobol LUT.\n', ...
            nGrid, cfg.lut.maxGridCandidates);
    end
end

skip0 = cfg.lut.sobolSkipOffset;
Grand = local_candidate_cloud(lb, ub, cfg.lut.nGlobal, skip0);
Gglobal = [Gprior; Ggrid; Grand];
sourceGlobal = [ones(size(Gprior,1),1,'uint8'); ...
    4*ones(size(Ggrid,1),1,'uint8'); 2*ones(size(Grand,1),1,'uint8')];

hitLookup = local_build_peak_hit_lookup(ref, cfg);
scorerCoarse = local_prepare_scorer(ptsCoarse, ref, det, cfg, computeInfo, hitLookup);
scorerFull = local_prepare_scorer(ptsFull, ref, det, cfg, computeInfo, hitLookup);

if strcmp(computeInfo.activeDevice, 'auto')
    nBench = min(size(Gglobal,1), cfg.compute.autoBenchmarkCandidates);
    computeInfo = local_auto_select_device(computeInfo, Gglobal(1:nBench,:), scorerCoarse, cfg);
end

fprintf('\nGlobal LUT: evaluating %d candidates with %d coarse pixels on %s.\n', ...
    size(Gglobal,1), numel(ptsCoarse.x), computeInfo.activeDevice);
[scoresGlobalApprox, computeInfo.activeDevice] = local_score_candidate_matrix( ...
    Gglobal, scorerCoarse, cfg, computeInfo.activeDevice, cfg.lut.progressEvery);

[~, orderApprox] = sort(scoresGlobalApprox, 'descend');
nRescoreGlobal = min(cfg.lut.globalFullRescoreTopK, numel(orderApprox));
idxGlobalExact = orderApprox(1:nRescoreGlobal);
fprintf('Exact full-pixel rescoring of top %d global candidates on CPU.\n', nRescoreGlobal);
scoresGlobalExact = local_score_candidate_matrix_exact_cpu( ...
    Gglobal(idxGlobalExact,:), scorerFull, cfg, cfg.lut.progressEvery);

[~, orderGlobalExact] = sort(scoresGlobalExact, 'descend');
idxSeedPool = idxGlobalExact(orderGlobalExact);

[halfWidthLevels, nSeedLevels, nPerSeedLevels, nRescoreLevels, minDistLevels] = ...
    local_lut_level_settings(cfg);
nLevels = size(halfWidthLevels,1);

nSeed0 = nSeedLevels(1);
idxSeeds = local_select_distinct_indices(Gglobal, idxSeedPool, nSeed0, ...
    halfWidthLevels(1,:), minDistLevels(1));
seedP = Gglobal(idxSeeds,:);

allP = Gglobal;
allScores = scoresGlobalApprox;
allSource = sourceGlobal;

Pexact = Gglobal(idxGlobalExact,:);
Sexact = scoresGlobalExact;
SexactSource = sourceGlobal(idxGlobalExact);

localPByLevel = cell(nLevels,1);
localApproxScoresByLevel = cell(nLevels,1);
localExactPByLevel = cell(nLevels,1);
localExactScoresByLevel = cell(nLevels,1);
seedPByLevel = cell(nLevels+1,1);
seedPByLevel{1} = seedP;

for iLevel = 1:nLevels
    hw = halfWidthLevels(iLevel,:);
    nPerSeed = nPerSeedLevels(iLevel);
    nSeedsHere = min(nSeedLevels(iLevel), size(seedP,1));
    seedP = seedP(1:nSeedsHere,:);

    Glocal = zeros(0,5);
    for k = 1:nSeedsHere
        p0 = seedP(k,:);
        lbLoc = max(lb, p0 - hw);
        ubLoc = min(ub, p0 + hw);
        skip = skip0 + 100000*iLevel + k*(nPerSeed + 4096);
        Glocal = [Glocal; local_candidate_cloud(lbLoc, ubLoc, nPerSeed, skip)];  
    end

    if isempty(Glocal)
        break
    end

    fprintf('Local LUT level %d/%d: %d candidates around %d seeds, half-width [%s], on %s.\n', ...
        iLevel, nLevels, size(Glocal,1), nSeedsHere, num2str(hw,' %.3g'), computeInfo.activeDevice);

    [scoresLocalApprox, computeInfo.activeDevice] = local_score_candidate_matrix( ...
        Glocal, scorerFull, cfg, computeInfo.activeDevice, cfg.lut.progressEvery);

    [~, orderLocalApprox] = sort(scoresLocalApprox, 'descend');
    nRescoreLocal = min(nRescoreLevels(iLevel), numel(orderLocalApprox));
    idxLocalExact = orderLocalApprox(1:nRescoreLocal);
    fprintf('Exact full-pixel rescoring of top %d candidates from local level %d.\n', ...
        nRescoreLocal, iLevel);
    scoresLocalExact = local_score_candidate_matrix_exact_cpu( ...
        Glocal(idxLocalExact,:), scorerFull, cfg, cfg.lut.progressEvery);

    PlevelExact = Glocal(idxLocalExact,:);
    [~, orderLevelExact] = sort(scoresLocalExact, 'descend');

    allP = [allP; Glocal];  
    allScores = [allScores; scoresLocalApprox];  
    allSource = [allSource; repmat(uint8(4+iLevel), size(Glocal,1), 1)];  

    Pexact = [Pexact; PlevelExact];  
    Sexact = [Sexact; scoresLocalExact];  
    SexactSource = [SexactSource; repmat(uint8(4+iLevel), numel(scoresLocalExact), 1)];  

    localPByLevel{iLevel} = Glocal;
    localApproxScoresByLevel{iLevel} = scoresLocalApprox;
    localExactPByLevel{iLevel} = PlevelExact;
    localExactScoresByLevel{iLevel} = scoresLocalExact;

    if iLevel < nLevels
        nSeedNext = nSeedLevels(iLevel+1);
        idxNext = local_select_distinct_indices(PlevelExact, orderLevelExact, nSeedNext, ...
            halfWidthLevels(iLevel+1,:), minDistLevels(iLevel+1));
        seedP = PlevelExact(idxNext,:);
        seedPByLevel{iLevel+1} = seedP;
    end
end

[Pexact, ia] = unique(round(Pexact, 10), 'rows', 'stable');
Sexact = Sexact(ia);
SexactSource = SexactSource(ia);
[sSortedExact, orderExact] = sort(Sexact, 'descend');
topK = min(cfg.lut.topK, numel(orderExact));

coarse = struct();
coarse.allP = allP;
coarse.allScores = allScores;
coarse.allSourceCode = allSource;
coarse.sortedScores = sort(allScores, 'descend');
coarse.exactP = Pexact;
coarse.exactScores = Sexact;
coarse.exactSourceCode = SexactSource;
coarse.exactOrder = orderExact;
coarse.sortedExactScores = sSortedExact;
coarse.topP = Pexact(orderExact(1:topK),:);
coarse.topScores = sSortedExact(1:topK);
coarse.lb = lb;
coarse.ub = ub;
coarse.seedPByLevel = seedPByLevel;
coarse.localPByLevel = localPByLevel;
coarse.localApproxScoresByLevel = localApproxScoresByLevel;
coarse.localExactPByLevel = localExactPByLevel;
coarse.localExactScoresByLevel = localExactScoresByLevel;
coarse.nGlobalCandidates = size(Gglobal,1);
coarse.nLocalCandidates = size(allP,1) - size(Gglobal,1);
coarse.computeDevice = computeInfo.activeDevice;

fprintf('\nTop exact LUT candidates:\n');
for k = 1:min(10, topK)
    p = coarse.topP(k,:);
    fprintf('  %2d score %.6f  p=[%+.2f %+.2f %+.2f  %+.2f %+.2f]\n', ...
        k, coarse.topScores(k), p);
end
end

% ======================================================================
function [halfWidths, nSeeds, nPerSeed, nRescore, minDistance] = local_lut_level_settings(cfg)
halfWidths = double(cfg.lut.localHalfWidth);
if isvector(halfWidths)
    if numel(halfWidths) ~= 5
        error('cfg.lut.localHalfWidth must have 5 columns.');
    end
    halfWidths = reshape(halfWidths,1,5);
elseif size(halfWidths,2) ~= 5
    error('cfg.lut.localHalfWidth must have 5 columns.');
end
nLevels = size(halfWidths,1);
nSeeds = local_repeat_level_value(cfg.lut.nLocalSeeds, nLevels, 'cfg.lut.nLocalSeeds');
nPerSeed = local_repeat_level_value(cfg.lut.nLocalPerSeed, nLevels, 'cfg.lut.nLocalPerSeed');
nRescore = local_repeat_level_value(cfg.lut.localFullRescoreTopK, nLevels, ...
    'cfg.lut.localFullRescoreTopK');
minDistance = local_repeat_level_value(cfg.lut.seedMinDistance, nLevels, ...
    'cfg.lut.seedMinDistance');

nSeeds = max(1, round(nSeeds));
nPerSeed = max(1, round(nPerSeed));
nRescore = max(1, round(nRescore));
minDistance = max(0, minDistance);
end

% ======================================================================
function v = local_repeat_level_value(vIn, nLevels, name)
v = double(vIn(:).');
if isscalar(v)
    v = repmat(v,1,nLevels);
elseif numel(v) ~= nLevels
    error('%s must be scalar or have one value per local LUT level.', name);
end
end

% ======================================================================
function lookup = local_build_peak_hit_lookup(ref, cfg)
dq = cfg.compute.hitLookupDQ_Ainv;
if ~isscalar(dq) || ~isfinite(dq) || dq <= 0
    error('cfg.compute.hitLookupDQ_Ainv must be a positive scalar.');
end
qMax = max([ref.qGrid_Ainv(:); ref.qPeaks_Ainv(:)]) + 5*cfg.lut.peakSigmaQ_Ainv;
qGrid = (0:dq:qMax).';
[dmin, ~] = local_nearest_q(qGrid, ref.qPeaks_Ainv, ref.qMidEdges_Ainv);
hit = exp(-0.5 * (dmin / cfg.lut.peakSigmaQ_Ainv).^2);
lookup = struct('dq', dq, 'invDQ', 1/dq, 'qMax', qGrid(end), 'hitDouble', hit(:));
end

% ======================================================================
function scorer = local_prepare_scorer(pts, ref, det, cfg, computeInfo, lookup)
w = double(pts.w(:));
w = w / max(sum(w), eps);
scorer = struct();
scorer.x = double(pts.x(:).');
scorer.y = double(pts.y(:).');
scorer.w = w;
scorer.pixel_mm = det.pixel_mm;
scorer.qScale = 4*pi/ref.lambda_A;
scorer.qRef = ref.qPeaks_Ainv(:);
scorer.qMidEdges = ref.qMidEdges_Ainv(:);
scorer.lookup = lookup;
scorer.nPoints = numel(pts.x);
scorer.hasGPU = false;

if computeInfo.gpuAvailable && ~strcmp(computeInfo.requestedDevice, 'cpu')
    try
        switch lower(char(cfg.compute.gpuPrecision))
            case 'double'
                castFcn = @double;
            otherwise
                castFcn = @single;
        end
        scorer.gx = gpuArray(castFcn(scorer.x));
        scorer.gy = gpuArray(castFcn(scorer.y));
        scorer.gw = gpuArray(castFcn(scorer.w));
        scorer.gHit = gpuArray(castFcn(lookup.hitDouble));
        scorer.gPixelMM = castFcn(scorer.pixel_mm);
        scorer.gQScale = castFcn(scorer.qScale);
        scorer.gInvDQ = castFcn(lookup.invDQ);
        scorer.gNLut = numel(lookup.hitDouble);
        scorer.hasGPU = true;
    catch ME
        warning('GPU scorer initialization failed; CPU will be used: %s', ME.message);
        scorer.hasGPU = false;
    end
end
end

% ======================================================================
function [scores, activeDevice] = local_score_candidate_matrix(G, scorer, cfg, activeDevice, progressEvery)
activeDevice = lower(char(activeDevice));
if strcmp(activeDevice, 'gpu')
    try
        scores = local_score_candidate_matrix_gpu(G, scorer, cfg, progressEvery);
        return
    catch ME
        warning('GPU LUT scoring failed; falling back to CPU: %s', ME.message);
        activeDevice = 'cpu';
    end
end
scores = local_score_candidate_matrix_cpu(G, scorer, cfg, progressEvery);
end

% ======================================================================
function scores = local_score_candidate_matrix_cpu(G, scorer, cfg, progressEvery)
n = size(G,1);
blockSize = max(1, round(cfg.compute.cpuBlockCandidates));
nBlocks = ceil(n/blockSize);
useParallel = local_use_parallel_cpu(cfg.compute.cpuParallel);
scorerCPU = local_strip_gpu_fields(scorer);

if useParallel && nBlocks > 1
    blockScores = cell(nBlocks,1);
    parfor ib = 1:nBlocks
        i1 = (ib-1)*blockSize + 1;
        i2 = min(n, ib*blockSize);
        blockScores{ib} = local_score_candidate_block_cpu(G(i1:i2,:), scorerCPU, cfg);
    end
    scores = vertcat(blockScores{:});
    if progressEvery > 0
        fprintf('  evaluated %6d / %6d, current best %.6f\n', n, n, max(scores));
    end
else
    scores = -Inf(n,1);
    nextReport = progressEvery;
    for ib = 1:nBlocks
        i1 = (ib-1)*blockSize + 1;
        i2 = min(n, ib*blockSize);
        scores(i1:i2) = local_score_candidate_block_cpu(G(i1:i2,:), scorerCPU, cfg);
        if progressEvery > 0 && (i2 >= nextReport || i2 == n)
            fprintf('  evaluated %6d / %6d, current best %.6f\n', i2, n, max(scores(1:i2)));
            nextReport = nextReport + progressEvery;
        end
    end
end
end

% ======================================================================
function scorerCPU = local_strip_gpu_fields(scorer)
scorerCPU = scorer;
names = fieldnames(scorerCPU);
remove = startsWith(names, 'g') | strcmp(names, 'hasGPU');
if any(remove)
    scorerCPU = rmfield(scorerCPU, names(remove));
end
end

% ======================================================================
function scores = local_score_candidate_block_cpu(Gb, scorer, cfg)
q = local_q_old_candidate_block(Gb, scorer.x, scorer.y, scorer.pixel_mm, scorer.qScale);
idx = round(q * scorer.lookup.invDQ) + 1;
inRange = idx >= 1 & idx <= numel(scorer.lookup.hitDouble);
idx = min(max(idx, 1), numel(scorer.lookup.hitDouble));
hit = scorer.lookup.hitDouble(idx);
hit(~inRange) = 0;
scores = hit * scorer.w;
scores = local_apply_tilt_prior(scores, Gb, cfg);
end

% ======================================================================
function scores = local_score_candidate_matrix_gpu(G, scorer, cfg, progressEvery)
if ~scorer.hasGPU
    error('GPU scorer was not initialized.');
end
n = size(G,1);
blockSize = max(1, round(cfg.compute.gpuBlockCandidates));
nBlocks = ceil(n/blockSize);
scores = -Inf(n,1);
nextReport = progressEvery;
gdev = gpuDevice;

for ib = 1:nBlocks
    i1 = (ib-1)*blockSize + 1;
    i2 = min(n, ib*blockSize);
    switch lower(char(cfg.compute.gpuPrecision))
        case 'double'
            gG = gpuArray(double(G(i1:i2,:)));
        otherwise
            gG = gpuArray(single(G(i1:i2,:)));
    end
    q = local_q_old_candidate_block_gpu(gG, scorer);
    idx = round(q .* scorer.gInvDQ) + 1;
    inRange = idx >= 1 & idx <= scorer.gNLut;
    idx = min(max(idx, 1), scorer.gNLut);
    hit = scorer.gHit(idx) .* inRange;
    s = hit * scorer.gw;
    s = gather(s);
    scores(i1:i2) = local_apply_tilt_prior(double(s), G(i1:i2,:), cfg);
    if progressEvery > 0 && (i2 >= nextReport || i2 == n)
        wait(gdev);
        fprintf('  evaluated %6d / %6d, current best %.6f\n', i2, n, max(scores(1:i2)));
        nextReport = nextReport + progressEvery;
    end
end
wait(gdev);
end

% ======================================================================
function scores = local_score_candidate_matrix_exact_cpu(G, scorer, cfg, progressEvery)
n = size(G,1);
blockSize = max(1, min(round(cfg.compute.cpuBlockCandidates), 512));
nBlocks = ceil(n/blockSize);
scores = -Inf(n,1);
nextReport = progressEvery;
for ib = 1:nBlocks
    i1 = (ib-1)*blockSize + 1;
    i2 = min(n, ib*blockSize);
    q = local_q_old_candidate_block(G(i1:i2,:), scorer.x, scorer.y, scorer.pixel_mm, scorer.qScale);
    [dq, ~] = local_nearest_q(q, scorer.qRef, scorer.qMidEdges);
    hit = exp(-0.5 * (dq / cfg.lut.peakSigmaQ_Ainv).^2);
    s = hit * scorer.w;
    scores(i1:i2) = local_apply_tilt_prior(s, G(i1:i2,:), cfg);
    if progressEvery > 0 && (i2 >= nextReport || i2 == n)
        fprintf('  rescored  %6d / %6d, current best %.6f\n', i2, n, max(scores(1:i2)));
        nextReport = nextReport + progressEvery;
    end
end
end

% ======================================================================
function q = local_q_old_candidate_block(Gb, xRow, yRow, pixel_mm, qScale)
cx = Gb(:,1);
cy = Gb(:,2);
z = Gb(:,3);
tx = Gb(:,4);
ty = Gb(:,5);

u = (xRow - cx) * pixel_mm;
v = (yRow - cy) * pixel_mm;
a = -sind(ty).*u + (cosd(ty).*sind(tx)).*v;
r2 = u.^2 + v.^2 + z.^2 + 2*z.*a;
cosAlpha = (z + a) ./ sqrt(max(r2, realmin('double')));
cosAlpha = max(min(cosAlpha, 1), -1);
q = qScale .* sqrt(max(0, 0.5*(1-cosAlpha)));
end

% ======================================================================
function q = local_q_old_candidate_block_gpu(Gb, scorer)
cx = Gb(:,1);
cy = Gb(:,2);
z = Gb(:,3);
tx = Gb(:,4);
ty = Gb(:,5);
u = (scorer.gx - cx) .* scorer.gPixelMM;
v = (scorer.gy - cy) .* scorer.gPixelMM;
a = -sind(ty).*u + (cosd(ty).*sind(tx)).*v;
r2 = u.^2 + v.^2 + z.^2 + 2*z.*a;
cosAlpha = (z + a) ./ sqrt(r2);
cosAlpha = max(min(cosAlpha, 1), -1);
q = scorer.gQScale .* sqrt(max(0, 0.5*(1-cosAlpha)));
end

% ======================================================================
function scores = local_apply_tilt_prior(scores, G, cfg)
scores = scores(:) - local_tilt_prior_penalty(G, cfg);
end

% ======================================================================
function penalty = local_tilt_prior_penalty(P, cfg)
P = double(P);
if isempty(P)
    penalty = zeros(0,1);
    return
end
if cfg.prior.tiltPenaltyWeight <= 0 || strcmpi(cfg.prior.penaltyType,'none')
    penalty = zeros(size(P,1),1);
    return
end

tiltMag = hypot(P(:,4), P(:,5));
d = (tiltMag - cfg.prior.tiltDeg) / max(cfg.prior.tiltSigma_deg, eps);
switch lower(char(cfg.prior.penaltyType))
    case 'soft'
        penalty = cfg.prior.tiltPenaltyWeight * asinh(d).^2;
    case 'quadratic'
        penalty = cfg.prior.tiltPenaltyWeight * d.^2;
    case 'none'
        penalty = zeros(size(d));
    otherwise
        error('cfg.prior.penaltyType must be soft, quadratic, or none.');
end
if isfinite(cfg.prior.maxPenalty)
    penalty = min(penalty, max(cfg.prior.maxPenalty,0));
end
end

% ======================================================================
function tf = local_use_parallel_cpu(flag)
tf = false;
if islogical(flag)
    if ~flag || exist('gcp','file') ~= 2
        return
    end
    try
        if isempty(gcp('nocreate'))
            parpool('threads');
        end
        tf = true;
    catch
        tf = false;
    end
elseif ischar(flag) || isstring(flag)
    if strcmpi(flag, 'auto') && exist('gcp','file') == 2
        try
            tf = ~isempty(gcp('nocreate'));
        catch
            tf = false;
        end
    elseif strcmpi(flag, 'true') || strcmpi(flag, 'on')
        tf = local_use_parallel_cpu(true);
    end
end
end

% ======================================================================
function info = local_auto_select_device(info, Gsample, scorer, cfg)
if ~info.gpuAvailable || ~scorer.hasGPU || ~cfg.compute.autoBenchmark
    if info.gpuAvailable && scorer.hasGPU
        info.activeDevice = 'gpu';
    else
        info.activeDevice = 'cpu';
    end
    return
end

% Warm the GPU outside the timed region.
nWarm = min(128, size(Gsample,1));
local_score_candidate_matrix_gpu(Gsample(1:nWarm,:), scorer, cfg, 0);
gdev = gpuDevice;
wait(gdev);

t0 = tic;
local_score_candidate_matrix_cpu(Gsample, scorer, cfg, 0);
tCpu = toc(t0);

t0 = tic;
local_score_candidate_matrix_gpu(Gsample, scorer, cfg, 0);
wait(gdev);
tGpu = toc(t0);

info.benchmarkCPU_s = tCpu;
info.benchmarkGPU_s = tGpu;
if tGpu * cfg.compute.autoGpuMinSpeedup < tCpu
    info.activeDevice = 'gpu';
else
    info.activeDevice = 'cpu';
end
fprintf('Auto device benchmark on %d candidates: CPU %.3f s, GPU %.3f s -> %s.\n', ...
    size(Gsample,1), tCpu, tGpu, info.activeDevice);
end

% ======================================================================
function idxSelected = local_select_distinct_indices(P, orderedIdx, nWanted, scale, minDistance)
idxSelected = zeros(0,1);
scale = max(abs(scale(:).'), eps);
for j = 1:numel(orderedIdx)
    idx = orderedIdx(j);
    if isempty(idxSelected)
        idxSelected = idx;
    else
        d = sqrt(sum(((P(idx,:) - P(idxSelected,:)) ./ scale).^2, 2));
        if all(d >= minDistance)
            idxSelected(end+1,1) = idx;  
        end
    end
    if numel(idxSelected) >= nWanted
        break
    end
end
if isempty(idxSelected) && ~isempty(orderedIdx)
    idxSelected = orderedIdx(1);
end
end

function fit = local_refine_from_lut(coarse, pts, ref, det, cfg)
Kmax = min(cfg.opt.topKToRefine, size(coarse.topP,1));
ordered = (1:size(coarse.topP,1)).';
[halfWidths, ~, ~, ~, ~] = local_lut_level_settings(cfg);
refineScale = halfWidths(end,:);
idxUse = local_select_distinct_indices(coarse.topP, ordered, Kmax, ...
    refineScale, cfg.opt.refineMinDistance);
K = numel(idxUse);

allP = nan(K,5);
allScore = -Inf(K,1);
allPeakScore = -Inf(K,1);
allPriorPenalty = Inf(K,1);
allAssign = cell(K,1);
allInfo = cell(K,1);

fprintf('\nRefining %d distinct exact-LUT candidates with robust q residuals.\n', K);
for k = 1:K
    p = coarse.topP(idxUse(k),:);
    infoCycles = cell(cfg.opt.reassignCycles,1);
    assignment = local_assign_points(p, pts, ref, cfg);

    for cyc = 1:cfg.opt.reassignCycles
        if nnz(assignment.keep) < cfg.assign.minAssignedPixels
            warning('Candidate %d cycle %d has only %d assigned pixels.', k, cyc, nnz(assignment.keep));
            break
        end

        pOld = p;
        qAssignedOld = assignment.qAssigned;
        resFcn = @(pp) local_old_residual_vector(pp, pts, assignment, det, cfg);
        [p, fval, info] = local_bounded_least_squares( ...
            resFcn, p, coarse.lb, coarse.ub, cfg.opt);  
        p = min(max(p, coarse.lb), coarse.ub);
        assignmentNew = local_assign_points(p, pts, ref, cfg);

        paramStep = max(abs((p-pOld) ./ max(coarse.ub-coarse.lb, eps)));
        assignmentSame = isequaln(assignmentNew.qAssigned, qAssignedOld);
        info.paramStepNormalized = paramStep;
        info.assignmentUnchanged = assignmentSame;
        infoCycles{cyc} = info;
        assignment = assignmentNew;

        if assignmentSame && paramStep < cfg.opt.reassignParamTol
            break
        end
    end

    [scoreBeforeWRMS, peakScore, priorPenalty] = ...
        local_peak_hit_score_exact(p, pts, ref, det, cfg);
    score = scoreBeforeWRMS - 0.15 * assignment.wrms_Ainv;
    allP(k,:) = p;
    allScore(k) = score;
    allPeakScore(k) = peakScore;
    allPriorPenalty(k) = priorPenalty;
    allAssign{k} = assignment;
    allInfo{k} = infoCycles;
    fprintf('  refine %2d score %.6f  peak %.6f  wrms %.5f  n %4d  p=[%+.3f %+.3f %+.3f  %+.3f %+.3f]\n', ...
        k, score, peakScore, assignment.wrms_Ainv, nnz(assignment.keep), p);
end

[~, ib] = max(allScore);
fit = struct();
fit.bestP = allP(ib,:);
fit.bestScore = allScore(ib);
fit.bestPeakScore = allPeakScore(ib);
fit.bestPriorPenalty = allPriorPenalty(ib);
fit.bestAssignment = allAssign{ib};
fit.allP = allP;
fit.allScores = allScore;
fit.allPeakScores = allPeakScore;
fit.allPriorPenalties = allPriorPenalty;
fit.allAssignments = allAssign;
fit.allInfo = allInfo;
fit.coarseIndicesRefined = idxUse;
end

% ======================================================================
function assignment = local_assign_points(p, pts, ref, cfg)
q = local_q_old_points_fast(p, pts.x, pts.y, cfg.detector.pixel_um/1000, ref.lambda_A);
[dq, idx] = local_nearest_q(q, ref.qPeaks_Ainv, ref.qMidEdges_Ainv);
keep = isfinite(q) & isfinite(dq) & dq <= cfg.assign.maxDQ_Ainv;

if nnz(keep) > cfg.assign.maxRefinePixels
    keepIdx = find(keep);
    [~, ord] = sort(dq(keepIdx), 'ascend');
    take = keepIdx(ord(1:cfg.assign.maxRefinePixels));
    keep = false(size(keep));
    keep(take) = true;
end
qAssigned = nan(size(q));
qAssigned(keep) = ref.qPeaks_Ainv(idx(keep));
res = q(keep) - qAssigned(keep);
w = pts.w(keep);
if isempty(res)
    wrms = Inf;
else
    wrms = sqrt(sum(w(:) .* res(:).^2) / max(sum(w), eps));
end

assignment = struct();
assignment.qCalc = q;
assignment.dqNearest = dq;
assignment.idxNearest = idx;
assignment.keep = keep;
assignment.qAssigned = qAssigned;
assignment.wrms_Ainv = wrms;
assignment.nAssigned = nnz(keep);
assignment.usedPeakQ = unique(qAssigned(keep));
end

% ======================================================================
function r = local_old_residual_vector(p, pts, assignment, det, cfg)
id = assignment.keep;
q = local_q_old_points_fast(p, pts.x(id), pts.y(id), det.pixel_mm, ...
    12.398419843320026 / cfg.beam.EkeV);
res = q(:) - assignment.qAssigned(id);
w = double(pts.w(id));
w = w(:) / max(sum(w), eps);
r = sqrt(w) .* asinh(res / cfg.opt.robustScaleQ_Ainv);

priorPenalty = local_tilt_prior_penalty(p, cfg);
r(end+1,1) = sqrt(max(priorPenalty,0));
end

% ======================================================================
function [score, rawPeakScore, priorPenalty] = local_peak_hit_score_exact(p, pts, ref, det, cfg)
q = local_q_old_points_fast(p, pts.x, pts.y, det.pixel_mm, ref.lambda_A);
if any(~isfinite(q))
    score = -Inf;
    rawPeakScore = -Inf;
    priorPenalty = Inf;
    return
end
[dq, ~] = local_nearest_q(q, ref.qPeaks_Ainv, ref.qMidEdges_Ainv);
hit = exp(-0.5 * (dq / cfg.lut.peakSigmaQ_Ainv).^2);
w = pts.w(:);
rawPeakScore = sum(w .* hit(:)) / max(sum(w), eps);
priorPenalty = local_tilt_prior_penalty(p, cfg);
score = rawPeakScore - priorPenalty;
end

% ======================================================================
function q = local_q_old_points_fast(p, xPix, yPix, pixel_mm, lambda_A)
% Exact simplified old APS geometry.  Only the z component of the rotated
% in-plane displacement is needed for the scattering angle.
cx = p(1);
cy = p(2);
z = p(3);
tx = p(4);
ty = p(5);

u = (double(xPix(:)) - cx) * pixel_mm;
v = (double(yPix(:)) - cy) * pixel_mm;
a = -sind(ty).*u + cosd(ty).*sind(tx).*v;
r2 = u.^2 + v.^2 + z.^2 + 2*z.*a;
cosAlpha = (z + a) ./ sqrt(max(r2, realmin('double')));
cosAlpha = max(min(cosAlpha, 1), -1);
q = (4*pi/lambda_A) .* sqrt(max(0, 0.5*(1-cosAlpha)));
end

% ======================================================================
function geom = local_build_geom_old(imSize, pixel_mm, p, EkeV, polMode)
nRows = imSize(1);
nCols = imSize(2);
[xPix, yPix] = meshgrid(1:nCols, 1:nRows);
cx = p(1); cy = p(2); z = p(3); tx_deg = p(4); ty_deg = p(5);
lambda_A = 12.398419843320026 / EkeV;

u = (xPix - cx) * pixel_mm;
v = (yPix - cy) * pixel_mm;
tx = deg2rad(tx_deg);
ty = deg2rad(ty_deg);
Rx = [1 0 0; 0 cos(tx) -sin(tx); 0 sin(tx) cos(tx)];
Ry = [cos(ty) 0 sin(ty); 0 1 0; -sin(ty) 0 cos(ty)];
Rmat = Ry * Rx;
eu = Rmat(:,1);
ev = Rmat(:,2);
n = Rmat(:,3);
X = u * eu(1) + v * ev(1);
Y = u * eu(2) + v * ev(2);
Z = z + u * eu(3) + v * ev(3);
Rpix = sqrt(X.^2 + Y.^2 + Z.^2);
sx = X ./ max(Rpix, eps);
sy = Y ./ max(Rpix, eps);
sz = Z ./ max(Rpix, eps);
cosAlpha = max(min(sz, 1), -1);
twoThetaMap_deg = acosd(cosAlpha);
qMap = (4*pi/lambda_A) .* sin(0.5 * acos(cosAlpha));
phiMap = atan2(Y, X);
ndotS = abs(n(1)*sx + n(2)*sy + n(3)*sz);
geom0 = z^2 / max(abs(n(3)), eps);
geometryCorr = (Rpix.^2 ./ max(ndotS, eps)) / max(geom0, eps);
switch lower(polMode)
    case 'none'
        polFactor = ones(size(qMap));
    case 'horizontal'
        polFactor = 1 - sx.^2;
    case 'vertical'
        polFactor = 1 - sy.^2;
    otherwise
        error('polMode must be horizontal, vertical, or none.');
end
polarizationCorr = 1 ./ max(polFactor, eps);
totalCorr = geometryCorr .* polarizationCorr;

geom = struct();
geom.qMap = qMap;
geom.phiMap = phiMap;
geom.twoThetaMap_deg = twoThetaMap_deg;
geom.geometryCorr = geometryCorr;
geom.polarizationCorr = polarizationCorr;
geom.totalCorr = totalCorr;
geom.center = [cx cy];
geom.z0_mm = z;
geom.z0_um = z * 1000;
geom.tx_deg = tx_deg;
geom.ty_deg = ty_deg;
geom.pixel_um = pixel_mm * 1000;
geom.EkeV = EkeV;
geom.polMode = polMode;
geom.Rmat = Rmat;
geom.normal = n;
end

% ======================================================================
function prof = local_profile_from_geom(Im, geom, qEdges, trimPct)
qEdges = qEdges(:).';
qCenters = qEdges(1:end-1) + 0.5*diff(qEdges);
Icorr = Im .* geom.totalCorr;
valid = isfinite(Im) & isfinite(Icorr) & isfinite(geom.qMap) & Im > 0 & ...
    geom.qMap >= qEdges(1) & geom.qMap < qEdges(end);
qVals = geom.qMap(valid);
rawVals = Im(valid);
corrVals = Icorr(valid);
[~,~,qbin] = histcounts(qVals, qEdges);
profileRaw = nan(numel(qCenters),1);
profileCorr = nan(numel(qCenters),1);
counts = zeros(numel(qCenters),1);
for k = 1:numel(qCenters)
    id = qbin == k;
    counts(k) = nnz(id);
    if counts(k) == 0
        continue
    end
    vr = rawVals(id);
    vr = vr(isfinite(vr) & vr > 0);
    vc = corrVals(id);
    vc = vc(isfinite(vc) & vc > 0);
    if ~isempty(vr)
        profileRaw(k) = local_trimmean(vr, trimPct);
    end
    if ~isempty(vc)
        profileCorr(k) = local_trimmean(vc, trimPct);
    end
end
prof = struct('qCenters', qCenters(:), 'profileRaw', profileRaw, 'profileCorr', profileCorr, 'counts', counts);
end

% ======================================================================
function pose6d = local_old_params_to_pose6d(p, imSize, pixel_mm)
% Convert beam-intercept old geometry to detector-center pose.
% R = Ry(ty)*Rx(tx), detector center is pixel [(nCols+1)/2, (nRows+1)/2].
nRows = imSize(1);
nCols = imSize(2);
cDet = [(nCols+1)/2, (nRows+1)/2];
cx = p(1); cy = p(2); z = p(3); tx_deg = p(4); ty_deg = p(5);
tx = deg2rad(tx_deg);
ty = deg2rad(ty_deg);
Rx = [1 0 0; 0 cos(tx) -sin(tx); 0 sin(tx) cos(tx)];
Ry = [cos(ty) 0 sin(ty); 0 1 0; -sin(ty) 0 cos(ty)];
R = Ry * Rx;
% Beam point has lab [0 0 z] and local detector coordinate of pixel [cx,cy].
uBeam = (cx - cDet(1)) * pixel_mm;
vBeam = (cy - cDet(2)) * pixel_mm;
T = [0;0;z] - R * [uBeam; vBeam; 0];
pose6d = [T(:).' tx_deg ty_deg 0];
end

% ======================================================================
function [lb, ub] = local_global_bounds(cfg)
lb = [cfg.bounds.beamX_px(1), cfg.bounds.beamY_px(1), cfg.bounds.zBeam_mm(1), ...
      cfg.bounds.thetaX_deg(1), cfg.bounds.thetaY_deg(1)];
ub = [cfg.bounds.beamX_px(2), cfg.bounds.beamY_px(2), cfg.bounds.zBeam_mm(2), ...
      cfg.bounds.thetaX_deg(2), cfg.bounds.thetaY_deg(2)];

% The main-axis prior guides the anchor grid and the soft score.  It changes
% the hard user bounds only when explicitly requested.
axisName = lower(char(cfg.prior.mainTiltAxis));
if ~any(strcmp(axisName, {'x','y','both'}))
    error('cfg.prior.mainTiltAxis must be x, y, or both.');
end
if cfg.prior.constrainBounds
    if strcmp(axisName, 'x')
        mag = cfg.prior.tiltDeg;
        if cfg.prior.tiltSign > 0
            lb(4) = max(lb(4), mag - 2*cfg.prior.tiltSigma_deg);
            ub(4) = min(ub(4), mag + 2*cfg.prior.tiltSigma_deg);
        elseif cfg.prior.tiltSign < 0
            lb(4) = max(lb(4), -mag - 2*cfg.prior.tiltSigma_deg);
            ub(4) = min(ub(4), -mag + 2*cfg.prior.tiltSigma_deg);
        end
    elseif strcmp(axisName, 'y')
        mag = cfg.prior.tiltDeg;
        if cfg.prior.tiltSign > 0
            lb(5) = max(lb(5), mag - 2*cfg.prior.tiltSigma_deg);
            ub(5) = min(ub(5), mag + 2*cfg.prior.tiltSigma_deg);
        elseif cfg.prior.tiltSign < 0
            lb(5) = max(lb(5), -mag - 2*cfg.prior.tiltSigma_deg);
            ub(5) = min(ub(5), -mag + 2*cfg.prior.tiltSigma_deg);
        end
    end
end
if any(lb >= ub)
    error('Invalid parameter bounds after applying prior. lb = [%s], ub = [%s]', num2str(lb), num2str(ub));
end
end

% ======================================================================

% ======================================================================
function G = local_prior_grid_candidates(lb, ub, cfg)
xVec = linspace(lb(1), ub(1), cfg.lut.priorGridXCount);
yVec = linspace(lb(2), ub(2), cfg.lut.priorGridYCount);
zVec = lb(3):cfg.lut.priorGridZStep_mm:ub(3);
if isempty(zVec) || abs(zVec(end)-ub(3)) > 10*eps
    zVec = unique([zVec ub(3)]);
end

axisName = lower(char(cfg.prior.mainTiltAxis));
main0 = cfg.prior.tiltSign * cfg.prior.tiltDeg;
mainOffsets = cfg.lut.priorMainTiltOffsets_deg;
minorOffsets = cfg.lut.priorMinorTiltOffsets_deg;
if strcmp(axisName, 'x')
    txVec = main0 + mainOffsets;
    tyVec = minorOffsets;
elseif strcmp(axisName, 'y')
    txVec = minorOffsets;
    tyVec = main0 + mainOffsets;
else
    txVec = lb(4):4:ub(4);
    tyVec = lb(5):4:ub(5);
end
if cfg.prior.tiltSign == 0
    if strcmp(axisName, 'x')
        txVec = unique([cfg.prior.tiltDeg + mainOffsets, -cfg.prior.tiltDeg + mainOffsets]);
    elseif strcmp(axisName, 'y')
        tyVec = unique([cfg.prior.tiltDeg + mainOffsets, -cfg.prior.tiltDeg + mainOffsets]);
    end
end

txVec = txVec(txVec >= lb(4) & txVec <= ub(4));
tyVec = tyVec(tyVec >= lb(5) & tyVec <= ub(5));
if isempty(txVec), txVec = (lb(4)+ub(4))/2; end
if isempty(tyVec), tyVec = (lb(5)+ub(5))/2; end
[X,Y,Z,TX,TY] = ndgrid(xVec, yVec, zVec, txVec, tyVec);
G = [X(:) Y(:) Z(:) TX(:) TY(:)];

if cfg.lut.includeBroadAnchorGrid
    xb = linspace(lb(1), ub(1), max(2,round(cfg.lut.broadAnchorXCount)));
    yb = linspace(lb(2), ub(2), max(2,round(cfg.lut.broadAnchorYCount)));
    zb = linspace(lb(3), ub(3), max(2,round(cfg.lut.broadAnchorZCount)));
    txb = local_anchor_vector(lb(4), ub(4), cfg.lut.broadAnchorTiltStep_deg);
    tyb = local_anchor_vector(lb(5), ub(5), cfg.lut.broadAnchorTiltStep_deg);
    [Xb,Yb,Zb,TXb,TYb] = ndgrid(xb,yb,zb,txb,tyb);
    G = [G; Xb(:) Yb(:) Zb(:) TXb(:) TYb(:)];  
end

G = min(max(G, lb), ub);
G = unique(round(G, 8), 'rows');
end

% ======================================================================
function v = local_anchor_vector(lo, hi, step)
step = max(abs(step), eps);
v = lo:step:hi;
v = unique([lo, v, hi, 0]);
v = v(v >= lo & v <= hi);
if isempty(v)
    v = 0.5*(lo+hi);
end
end

% ======================================================================
function G = local_candidate_cloud(lb, ub, n, skipOffset)
D = numel(lb);
if n <= 0
    G = zeros(0,D);
    return
end
if nargin < 4
    skipOffset = 0;
end
if exist('sobolset','file') == 2
    try
        s = sobolset(D, 'Skip', 1024 + skipOffset, 'Leap', 97);
        s = scramble(s, 'MatousekAffineOwen');
        U = net(s, n);
    catch
        U = rand(n, D);
    end
else
    U = rand(n, D);
end
G = lb + U .* (ub - lb);
end

% ======================================================================
function [dmin, idx] = local_nearest_q(q, qRef, midEdges)
% Exact nearest neighbor for a sorted one-dimensional peak list.
qRef = qRef(:);
if nargin < 3 || isempty(midEdges)
    midEdges = [-Inf; 0.5*(qRef(1:end-1) + qRef(2:end)); Inf];
end
idx = discretize(q, midEdges);
qNearest = qRef(idx(:));
qNearest = reshape(qNearest, size(q));
dmin = abs(q - qNearest);
end

% ======================================================================
function y = local_gaussian_fwhm(x, x0, fwhm)
sigma = fwhm / (2*sqrt(2*log(2)));
y = exp(-0.5*((x(:)-x0)/sigma).^2);
end

% ======================================================================
function [xBest, fBest, info] = local_bounded_least_squares(resFcn, x0, lb, ub, opt)
x0 = min(max(x0(:).', lb), ub);
span = ub - lb;
t0 = (x0 - lb) ./ span;
t0 = min(max(t0, 0), 1);
resT = @(t) resFcn(lb + t(:).' .* span);
solver = lower(char(opt.solver));
if strcmp(solver, 'auto')
    if exist('lsqnonlin','file') == 2
        solver = 'lsqnonlin';
    elseif exist('fmincon','file') == 2
        solver = 'fmincon';
    else
        solver = 'fminsearch';
    end
end

info = struct();
switch solver
    case 'lsqnonlin'
        if exist('lsqnonlin','file') ~= 2
            error('lsqnonlin was requested but is unavailable.');
        end
        options = optimoptions('lsqnonlin', 'Display','off', ...
            'Algorithm','trust-region-reflective', ...
            'MaxIterations', opt.maxIter, ...
            'MaxFunctionEvaluations', opt.maxFunEvals, ...
            'FiniteDifferenceType','forward', ...
            'StepTolerance', opt.stepTolerance, ...
            'FunctionTolerance', opt.functionTolerance, ...
            'OptimalityTolerance', opt.optimalityTolerance);
        [tBest, fBest, ~, exitflag, output] = lsqnonlin(resT, t0, zeros(size(t0)), ones(size(t0)), options);
        info.method = 'lsqnonlin';
    case 'fmincon'
        if exist('fmincon','file') ~= 2
            error('fmincon was requested but is unavailable.');
        end
        scalarObj = @(t) sum(resT(t).^2);
        options = optimoptions('fmincon', 'Display','off', 'Algorithm','interior-point', ...
            'MaxIterations', opt.maxIter, 'MaxFunctionEvaluations', opt.maxFunEvals, ...
            'StepTolerance', opt.stepTolerance, ...
            'FunctionTolerance', opt.functionTolerance, ...
            'OptimalityTolerance', opt.optimalityTolerance);
        [tBest, fBest, exitflag, output] = fmincon(scalarObj, t0, [], [], [], [], ...
            zeros(size(t0)), ones(size(t0)), [], options);
        info.method = 'fmincon';
    case {'fminsearch','nelder-mead'}
        scalarObj = @(t) sum(resT(t).^2);
        options = optimset('Display','off', 'MaxIter', opt.maxIter, ...
            'MaxFunEvals', opt.maxFunEvals, 'TolX', opt.stepTolerance, ...
            'TolFun', opt.functionTolerance);
        [tBest, fBest, exitflag, output] = local_fminsearchbnd( ...
            scalarObj, t0, zeros(size(t0)), ones(size(t0)), options);
        info.method = 'fminsearchbnd';
    otherwise
        error('cfg.opt.solver must be auto, lsqnonlin, fmincon, or fminsearch.');
end
xBest = lb + tBest(:).' .* span;
info.exitflag = exitflag;
info.output = output;
end

% ======================================================================
function [x, fval, exitflag, output] = local_fminsearchbnd(fun, x0, LB, UB, options)
% Lightweight bounded Nelder-Mead wrapper using sine transform.
x0 = x0(:).'; LB = LB(:).'; UB = UB(:).';
y0 = zeros(size(x0));
for i = 1:numel(x0)
    if isfinite(LB(i)) && isfinite(UB(i))
        t = 2*(x0(i)-LB(i))/(UB(i)-LB(i)) - 1;
        t = max(min(t, 0.999999), -0.999999);
        y0(i) = asin(t);
    else
        y0(i) = x0(i);
    end
end
wrap = @(y) fun(local_y_to_x(y, LB, UB));
[y, fval, exitflag, output] = fminsearch(wrap, y0, options);
x = local_y_to_x(y, LB, UB);
end

function x = local_y_to_x(y, LB, UB)
y = y(:).';
x = y;
for i = 1:numel(y)
    if isfinite(LB(i)) && isfinite(UB(i))
        x(i) = LB(i) + (UB(i)-LB(i)) * (sin(y(i)) + 1) / 2;
    elseif isfinite(LB(i))
        x(i) = LB(i) + exp(y(i));
    elseif isfinite(UB(i))
        x(i) = UB(i) - exp(y(i));
    end
end
end

% ======================================================================
function tripCfg = local_make_tripol_cfg(cfg, validMask)
tripCfg = cfg.tripol;
tripCfg.mask = validMask;
if ~isfield(tripCfg, 'varMap')
    tripCfg.varMap = [];
end
if ~isfield(tripCfg, 'phi0_deg')
    tripCfg.phi0_deg = 0;
end
end

% ======================================================================
function local_plot_signal_processing(out)
prep = out.prep;
pts = out.scorePoints;
p = out.params;
figure('Name','scan2 signal processing audit: old+LUT', 'Color','w', 'Position',[50 80 1500 850]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

nexttile
imagesc(prep.logImage); axis image ij; colorbar
hold on; plot(p(1), p(2), '-rx', 'MarkerSize',12, 'LineWidth',2)
title('log(1 + raw image), fitted beam intercept')
xlabel('pixel x'); ylabel('pixel y')

nexttile
imagesc(prep.background); axis image ij; colorbar
title(sprintf('moving-median background, window %d px', prep.backgroundWindowPx))
xlabel('pixel x'); ylabel('pixel y')

nexttile
imagesc(prep.enhanced); axis image ij; colorbar
title('ring-enhanced image = log image - background')
xlabel('pixel x'); ylabel('pixel y')

nexttile
vals = prep.enhanced(prep.validMask & isfinite(prep.enhanced));
histogram(vals, 120); set(gca,'YScale','log'); grid on; hold on
xline(prep.thresholdSigma, '--', 'sigma threshold');
xline(prep.thresholdPercentile, '--', 'percentile threshold');
xline(prep.threshold, 'k-', 'chosen threshold', 'LineWidth',1.5);
title('threshold audit: enhanced image, not raw-count threshold')
xlabel('enhanced intensity'); ylabel('pixel count')

nexttile
imagesc(prep.ringMask); axis image ij; colorbar
title(sprintf('ring-enhanced mask, %d pixels', nnz(prep.ringMask)))
xlabel('pixel x'); ylabel('pixel y')

nexttile
imagesc(prep.logImage); axis image ij; colorbar; hold on
plot(pts.x, pts.y, '.', 'LineStyle','none', 'MarkerSize', 4)
plot(p(1), p(2), '-rx', 'MarkerSize',12, 'LineWidth',2)
title(sprintf('score/refine pixels over raw image, %d pixels', numel(pts.x)))
xlabel('pixel x'); ylabel('pixel y')
end

% ======================================================================
function local_plot_pose_result(out)
geom = out.geom;
ref = out.ref;
prep = out.prep;
p = out.params;
fit = out.fit;
qMap = geom.qMap;
levels = ref.qPeaks_Ainv(ref.qPeaks_Ainv >= min(qMap(:)) & ref.qPeaks_Ainv <= max(qMap(:)));
if numel(levels) > out.cfg.output.nContourPeaks
    % Prefer the lower-Q peaks visible in this detector image.
    levels = levels(1:out.cfg.output.nContourPeaks);
end

figure('Name','scan2 old+LUT pose result', 'Color','w', 'Position',[40 60 1550 900]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

nexttile
imagesc(prep.logImage); axis image ij; colorbar; hold on
if ~isempty(levels)
    contour(qMap, levels, 'w', 'LineWidth', 0.9);
end
plot(p(1), p(2), '-rx', 'MarkerSize',12, 'LineWidth',2)
title('LaB6 Q contours from fitted old+LUT pose')
xlabel('pixel x'); ylabel('pixel y')

nexttile
imagesc(qMap); axis image ij; colorbar
title('fitted Q map (A^{-1})')
xlabel('pixel x'); ylabel('pixel y')

nexttile
imagesc(geom.twoThetaMap_deg); axis image ij; colorbar
title('fitted 2theta map (deg)')
xlabel('pixel x'); ylabel('pixel y')

nexttile
% Direct peak-position comparison.  The enhanced radial profile already has
% the slowly varying 2-D image background removed.  Here we remove only one
% constant residual floor, then determine one multiplicative reference scale
% from narrow windows around the expected LaB6 peaks.  No curved baseline,
% Q shift, or independent per-peak scaling is applied.
q = out.profileEnhanced.qCenters(:);
yObs = out.profileEnhanced.profileRaw(:);
nPerBin = out.profileEnhanced.counts(:);
yRefUnit = interp1(ref.qGrid_Ainv(:), ref.profile(:), q, 'linear', 0);

good = isfinite(q) & isfinite(yObs) & isfinite(yRefUnit) & nPerBin > 0;

% Remove only a constant floor.  This avoids the artificial curved baseline
% produced by a 1-D baseline model while leaving all peak shapes unchanged.
if any(good)
    obsFloor = local_prctile(yObs(good), 5);
else
    obsFloor = 0;
end
yObs = yObs - obsFloor;
yObs(yObs < 0) = 0;
yObs(~good) = NaN;

% Use all retained LaB6 peaks visible on this Q grid.
qPeaksPlot = ref.qPeaks_Ainv(:);
peakIPlot = ref.peakI(:);
keepPeak = isfinite(qPeaksPlot) & qPeaksPlot >= min(q) & qPeaksPlot <= max(q);
qPeaksPlot = qPeaksPlot(keepPeak);
peakIPlot = peakIPlot(keepPeak);

% Fit one nonnegative reference scale using only narrow peak windows, so the
% residual background between peaks does not determine the visual scaling.
peakHalfWidth_Ainv = 0.05;
scaleMask = false(size(q));
for k = 1:numel(qPeaksPlot)
    scaleMask = scaleMask | abs(q - qPeaksPlot(k)) <= peakHalfWidth_Ainv;
end
scaleMask = scaleMask & good & yRefUnit > 0.01*max(yRefUnit);

if nnz(scaleMask) >= 3
    scaleWeight = sqrt(max(nPerBin(scaleMask),1));
    refScale = sum(scaleWeight .* yObs(scaleMask) .* yRefUnit(scaleMask)) / ...
        max(sum(scaleWeight .* yRefUnit(scaleMask).^2), eps);
    refScale = max(refScale, 0);
else
    refScale = max(yObs(good)) / max(max(yRefUnit(good)), eps);
end
yRefScaled = refScale * yRefUnit;

% One common final normalization keeps both curves on the same plotted scale.
commonScale = max([yObs(good); yRefScaled(good); eps]);
yObsPlot = yObs / commonScale;
yRefPlot = yRefScaled / commonScale;

hRef = plot(q, yRefPlot, 'LineWidth',1.2); hold on
hObs = plot(q, yObsPlot, 'LineWidth',1.2);

% Short sticks show every expected LaB6 peak position without covering the
% curves.  Their modest height indicates simulated relative intensity.
hPeaks = gobjects(1);
if ~isempty(qPeaksPlot)
    if max(peakIPlot) > 0
        peakIPlot = peakIPlot / max(peakIPlot);
    else
        peakIPlot = ones(size(peakIPlot));
    end
    stickHeight = 0.025 + 0.055*sqrt(peakIPlot);
    for k = 1:numel(qPeaksPlot)
        hTmp = plot([qPeaksPlot(k) qPeaksPlot(k)], [0 stickHeight(k)], ':', ...
            'Color',[0.35 0.35 0.35], 'LineWidth',0.8);
        if k == 1
            hPeaks = hTmp;
        else
            set(hTmp, 'HandleVisibility','off');
        end
    end
end

xlim([min(q) max(q)]);
ylim([0 1.05]);
grid on
xlabel('Q (A^{-1})'); ylabel('common normalized intensity')
title('1D peak-position check: constant-offset removed, peak-scaled')
if isgraphics(hPeaks)
    legend([hRef hObs hPeaks], ...
        {'LaB6 reference, one fitted scale', ...
         'observed enhanced profile', ...
         'LaB6 peak positions'}, ...
        'Location','northeast');
else
    legend([hRef hObs], ...
        {'LaB6 reference, one fitted scale', ...
         'observed enhanced profile'}, ...
        'Location','northeast');
end

nexttile
s = out.coarse.sortedScores;
plot(s, '-', 'MarkerSize', 4); hold on
plot(numel(s)+(1:numel(fit.allScores)), fit.allScores, 'o', 'MarkerSize', 6)
grid on
xlabel('candidate rank'); ylabel('score')
title('candidate scores: LUT and refined candidates')
legend({'LUT sorted scores','refined top candidates'}, 'Location','best')

nexttile
axis off
normal = geom.normal(:).';
txt = sprintf(['p = [beamX beamY z tx ty]\n' ...
    '[%+.3f  %+.3f  %+.3f  %+.3f  %+.3f]\n\n' ...
    'pose6d = [Tx Ty Tz rx ry rz]\n' ...
    '[%+.3f  %+.3f  %+.3f  %+.3f  %+.3f  %+.3f]\n\n' ...
    'score = %.6f\npeak score = %.6f\nq WRMS = %.5f A^{-1}\n' ...
    'assigned pixels = %d\nsearch rounds = %d\nnormal = [%+.4f %+.4f %+.4f]'], ...
    out.params, out.pose6d, out.fit.bestScore, out.fit.bestPeakScore, ...
    out.fit.bestAssignment.wrms_Ainv, nnz(out.fit.bestAssignment.keep), ...
    numel(out.searchHistory), normal);
text(0.02, 0.95, txt, 'Units','normalized', 'VerticalAlignment','top', ...
    'FontName','FixedWidth', 'FontSize', 10)
title('fit summary')
end

% ======================================================================
function info = local_resolve_compute(opts)
requested = lower(char(opts.device));
if ~ismember(requested, {'auto','cpu','gpu'})
    error('cfg.compute.device must be auto, cpu, or gpu.');
end

gpuAvailable = false;
gpuName = '';
try
    if exist('canUseGPU','file') == 2
        gpuAvailable = canUseGPU;
    elseif exist('gpuDeviceCount','file') == 2
        gpuAvailable = gpuDeviceCount("available") > 0;
    end
    if gpuAvailable
        if opts.resetGPU
            reset(gpuDevice);
        end
        gd = gpuDevice;
        gpuName = gd.Name;
    end
catch
    gpuAvailable = false;
    gpuName = '';
end

if strcmp(requested, 'gpu') && ~gpuAvailable
    if opts.gpuRequired
        error('GPU was requested but no usable MATLAB GPU was found.');
    end
    warning('GPU was requested but unavailable; using CPU.');
    active = 'cpu';
elseif strcmp(requested, 'gpu')
    active = 'gpu';
elseif strcmp(requested, 'cpu')
    active = 'cpu';
elseif gpuAvailable
    active = 'auto';
else
    active = 'cpu';
end

info = struct();
info.requestedDevice = requested;
info.activeDevice = active;
info.gpuAvailable = gpuAvailable;
info.gpuName = gpuName;
info.benchmarkCPU_s = NaN;
info.benchmarkGPU_s = NaN;
end

% ======================================================================
function validation = local_validate_result(p, fit, opts)
validation = struct('enabled', false, 'passed', true, 'delta', [], 'tolerance', []);
if ~isfield(opts, 'referenceP') || isempty(opts.referenceP)
    return
end
refP = opts.referenceP(:).';
if numel(refP) ~= 5
    error('cfg.validation.referenceP must contain five values.');
end
tol = opts.tolerance(:).';
if isscalar(tol)
    tol = repmat(tol, 1, 5);
end
if numel(tol) ~= 5
    error('cfg.validation.tolerance must be scalar or contain five values.');
end
delta = p - refP;
passed = all(abs(delta) <= tol);
validation.enabled = true;
validation.passed = passed;
validation.referenceP = refP;
validation.delta = delta;
validation.tolerance = tol;
validation.score = fit.bestScore;
validation.wrms_Ainv = fit.bestAssignment.wrms_Ainv;

fprintf('Validation delta [beamX beamY z tx ty] = [%+.4g %+.4g %+.4g %+.4g %+.4g]\n', delta);
if ~passed
    msg = sprintf('Fast fit differs from the supplied reference beyond tolerance. Delta = [%s].', num2str(delta, ' %.5g'));
    if opts.warnOnly
        warning('%s', msg);
    else
        error('%s', msg);
    end
end
end

% ======================================================================
function local_save_output(outIn, opts)
out = outIn;   
if opts.saveCompact
    out.inputImage = [];
    out.imageUsed = [];
    if isfield(out, 'prep')
        out.prep.raw = [];
        out.prep.clipped = [];
    end
    if isfield(out, 'coarse')
        out.coarse.allP = [];
        out.coarse.allScores = [];
        out.coarse.allSourceCode = [];
    end
end
version = char(opts.saveVersion);
if isempty(version)
    version = '-v7.3';
end
save(opts.matFile, 'out', version);
end

% ======================================================================
function cfg = local_merge_struct(cfg, user)
if isempty(user)
    return
end
names = fieldnames(user);
for i = 1:numel(names)
    k = names{i};
    if isstruct(user.(k)) && isfield(cfg, k) && isstruct(cfg.(k))
        cfg.(k) = local_merge_struct(cfg.(k), user.(k));
    else
        cfg.(k) = user.(k);
    end
end
end

% ======================================================================
function f = local_resolve_file(primary, alternate)
f = '';
if exist(primary, 'file') == 2
    f = primary;
elseif nargin >= 2 && ~isempty(alternate) && exist(alternate, 'file') == 2
    f = alternate;
end
end

% ======================================================================
function p = local_prctile(x, pct)
x = x(isfinite(x));
if isempty(x)
    p = NaN;
    return
end
try
    p = prctile(x, pct);
catch
    x = sort(x(:));
    idx = 1 + (numel(x)-1) * pct/100;
    lo = floor(idx); hi = ceil(idx);
    if lo == hi
        p = x(lo);
    else
        p = x(lo) + (idx-lo) * (x(hi)-x(lo));
    end
end
end

% ======================================================================
function y = local_trimmean(x, pct)
x = x(isfinite(x));
if isempty(x)
    y = NaN;
    return
end
try
    y = trimmean(x, pct);
catch
    x = sort(x(:));
    n = numel(x);
    k = floor((pct/100) * n / 2);
    if 2*k >= n
        y = mean(x);
    else
        y = mean(x(k+1:n-k));
    end
end
end
