function out = cif2powder_CM(cifFile, EkeV, opts)
% This function reads a crystal structure from a CIF file and simulates
% its powder diffraction pattern for a chosen X-ray photon energy.
%
% It computes:
%   1) the discrete Bragg reflections ("stick pattern")
%   2) a broadened powder profile vs 2theta
%   3) a broadened powder profile vs Q
%
% The code is organized so that most physical corrections and profile
% choices can be turned on or off through the "opts" structure.
%
% Physical idea
% -------------
% A powder pattern is the orientational average of diffraction from many
% crystallites with random orientations.
%
% For each allowed reflection (h,k,l), the code:
%   - computes the d spacing from the unit cell metric
%   - computes the Bragg angle from the wavelength
%   - computes the structure factor F(hkl)
%   - forms an intensity proportional to abs(F)^2
%   - optionally applies standard diffraction corrections
%   - broadens the discrete lines into a continuous profile
%
% Peak positions are determined mainly by the lattice geometry.
% Peak intensities depend on the atomic scattering factors, occupancies,
% thermal motion, and the optional diffraction corrections.
%
%
% Inputs
% ------
% cifFile - Path to the CIF file containing the crystal structure.
%
% EkeV - X-ray photon energy in keV.
%
% opts - structure containing optional settings. Missing fields are filled
%        automatically with defaults. (see more later)
%
%
% Output
% ------
% out - structure containing the parsed crystal information, raw reflection
%        list, merged powder lines, and broadened profiles.
%
%   out.lambda_A - wavelength in Angstrom.
%   out.energy_keV - photon energy in keV.
%   out.cell - unit cell information read from the CIF.
%   out.atomsExpanded - atomic positions after expansion of the asymmetric
%                       unit by the CIF symmetry operations.
%   out.linesRaw - all enumerated reflections before powder-line merging.
%   out.linesMerged - powder lines merged by nearly equal Q.
%   out.tthGrid_deg - 2theta grid in degrees.
%   out.tthProfile -  simulated broadened intensity on the 2theta grid.
%   out.tthProfileNorm - same profile, normalized to max = 1.
%   out.qGrid_Ainv - q grid in A^-1.
%   out.qProfile - simulated broadened intensity on the Q grid.
%   out.qProfileNorm - same profile, normalized to max = 1.
%
% Coordinate and scattering conventions
% -------------------------------------
% Fractional coordinates from the CIF are used for the structure factor.
% Reciprocal-space convention:
%   Q = 4*pi*sin(theta)/lambda = 2*pi/d
% where:
%   theta   = Bragg angle
%   2theta  = diffraction angle
%   lambda  = wavelength
%   d       = lattice plane spacing
%
% -------------------------------------------------------------------------
% OPTIONS IN "opts"
% -------------------------------------------------------------------------
%
% The options are grouped into:
%
%   opts.twoThetaGrid_deg
%   opts.qGrid_Ainv
%   opts.corrections
%   opts.profile
%
%
% 1) Grid options
% -------------------------------------------------------------------------
%
% opts.twoThetaGrid_deg - vector for the output 2theta grid (degrees).
%   Only affects the sampled broadened profile, doesnt change Bragg peaks.
%
% opts.qGrid_Ainv - vector specifying the output Q grid in A^-1.
%                   only affect the sampling of the final profile.
%
% opts.mergeQtol_Ainv - tolerance in q used to merge reflections into powder lines.
%   Several reflections may land at the same or nearly the same powder
%   position. After computing all raw reflections, the code groups lines
%   whose Q values are equal within this tolerance.
%   If too small, reflections that should merge may remain split numerically.
%   If too large, nearby but distinct lines may merge incorrectly.
%
% opts.hklMax - optional manual control over the HKL search limits.
%   if left empty, the code builds a conservative rectangular search box in
%   HKL and then keeps only reflections satisfying the metric condition.
%   This parameter is mainly for debugging, speed control, or unusual cases.
%   Accepted forms: scalar    -> use same max for h, k, l
%                   3-vector  -> [hMax kMax lMax]
%
% opts.verbose - logical flag controlling warning and diagnostic messages.
%
%
% 2) Intensity correction options
% -------------------------------------------------------------------------
% opts.corrections.useOccupancy - if true, multiply each atomic contribution
%   by the site occupancy read from the CIF. Physical meaning: a partially
%   occupied site contributes less scattering than a fully occupied site.
%   If false, all occupancies are treated as 1.
%
% opts.corrections.useDebyeWaller - if true, use isotropic thermal motion
%      through the Debye-Waller factor. The code uses Biso from the CIF, or
%       converts Uiso to Biso when needed. The factor is exp(-B*s^2) where:
%       s = sin(theta)/lambda = Q/(4*pi)
%       capturing thermal motion reduces coherent Bragg intensity for
%       larger Q. If false, all B values are set to zero.
%
% opts.corrections.useLorentz - if true, apply the Lorentz correction:
%                                           1 / (sin(theta)^2 * cos(theta))
%   this changes relative intensities across angle.
%
% opts.corrections.usePolarization - if true, apply a standard unpolarized
% polarization factor: (1 + cos(2theta)^2) / 2 as the scattered intensity
% depends on polarization geometry.
%
% opts.corrections.useAnomalous - controls whether anomalous corrections are included.
%       if false   -> do not use anomalous corrections
%       if true    -> use anomalous corrections if anomalousFcn is supplied
%       if 'auto'  -> warn if the energy appears near an absorption edge
% as near an absorption edge, the atomic form factor becomes:
%                                              f = f0(Q) + fp(E) + i*fpp(E)
%
% opts.corrections.anomalousFcn
%   optional function handle returning [fp, fpp] for a given element and
%   energy:  [fp, fpp] = anomalousFcn(elemSymbol, EkeV)
%
% opts.corrections.nearEdgeWindow_keV - energy window used only for
%                        utomatic warning logic when useAnomalous = 'auto'.
%
%
% 3) Peak-profile options
% -------------------------------------------------------------------------
%
% opts.profile.widthModel - controls how the peak width depends on angle.
%       if 'fixed' all peaks have the same FWHM in 2theta
%          'caglioti' width depends on angle through the Caglioti equation
%
% opts.profile.shape - controls the functional form of the broadened line.
%       'gaussian' for symmetric Gaussian line shape
%       'pvoigt' for pseudo-Voigt (weighted sum of Gaussian & Lorentzian)
%
% opts.profile.fixedFWHM_deg - used only when widthModel = 'fixed'.
%   for the full width at half maximum (FWHM) of every peak in 2theta.
%
% opts.profile.U - broadening that grows like tan(theta)^2
% opts.profile.V - intermediate linear tan(theta) term
% opts.profile.W - constant baseline broadening
%      U/V/W Used only when widthModel = 'caglioti', these are the Caglioti
%      parameters in: FWHM(2theta)^2 = U*tan(theta)^2 + V*tan(theta) + W
%      where theta is the Bragg angle and FWHM is in degrees of 2theta.
%      These are empirical instrumental-profile parameters.
%
%
% opts.profile.useSize - if true, add crystallite-size broadening.
% opts.profile.size_A - effective crystallite size in Ang, when useSize=true.
%
% opts.profile.shapeK - Scherrer shape factor used in size broadening.
%                      adds a Scherrer-like contribution:
%                                     beta_size ~ K*lambda / (L*cos(theta))
%
% opts.profile.useStrain - if true, add microstrain broadening.
% opts.profile.microstrain - RMS microstrain used when useStrain = true.
% opts.profile.eta - Pseudo-Voigt mixing parameter, used when 'pvoigt' used.
%   eta = 0 pure Gaussian, eta = 1 pure Lorentzian
%
%
% Notes for new users
% -------------------
% 1) Peak positions vs intensities - the lattice controls peak positions.
%    The scattering factors and correction choices mainly control intensities.
%
% 2) Q pattern vs 2theta pattern - the Q positions of reflections are
%    determined by the crystal structure. The 2theta positions depend on
%    wavelength (hence photon energy EkeV).
%
% 3) Defaults are practical starting points - the default width and profile
%    values are not guaranteed to represent your experimental data.
%
% 4) If you care about a realistic simulation near an absorption edge
%    supply a real anomalous correction function for fp(E) and fpp(E).

% Example - copy, uncomment and run to see how it behaves:
% -------------------
% opts = struct;
% opts.twoThetaGrid_deg = 0.1:0.02:90;
% opts.profile.widthModel = 'caglioti';
% opts.profile.shape = 'pvoigt';
% opts.profile.U = 0.002;
% opts.profile.V = 0;
% opts.profile.W = 0.003;
% opts.profile.eta = 0.35;
% opts.corrections.useLorentz = true;
% opts.corrections.usePolarization = false;
%
% out = cif2powder_CM('1000057.cif', 20, opts); % for LaB6
%
% figure('Color','w');
% tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
%
% nexttile
% plot(out.tthGrid_deg, out.tthProfileNorm, 'k-', 'LineWidth', 1.3); hold on
% stem(out.linesMerged.twoTheta_deg, ...
%      out.linesMerged.I./max(out.linesMerged.I), ...
%      'r.', 'MarkerSize', 8)
% xlabel('2theta (deg)')
% ylabel('Normalized intensity')
% title('Powder diffraction profile vs 2theta')
% box on
% xlim([min(out.tthGrid_deg) max(out.tthGrid_deg)])
%
% nexttile
% plot(out.qGrid_Ainv, out.qProfileNorm, 'k-', 'LineWidth', 1.3); hold on
% stem(out.linesMerged.Q_Ainv, ...
%      out.linesMerged.I./max(out.linesMerged.I), ...
%      'r.', 'MarkerSize', 8)
% xlabel('Q (A^-1)')
% ylabel('Normalized intensity')
% title('Powder diffraction profile vs Q')
% box on
% xlim([min(out.qGrid_Ainv) max(out.qGrid_Ainv)])
%
%--------------------------------------------------------------------------%
%   Ver 1.3 (2022-03-12)
%   Adi Natan (natan@stanford.edu)
%--------------------------------------------------------------------------
%
% actual code starts here:


% Fill missing options with defaults.
if nargin < 3 || isempty(opts)
    opts = struct;
end
opts = set_default_opts_local(opts);

% Convert photon energy in keV to wavelength in Angstrom.
lambda_A = 12.398419843320026 / EkeV;

% Read the CIF and expand the asymmetric unit using symmetry.
cif = read_cif_basic_local(cifFile);
atoms = expand_asymmetric_unit_local(cif);

% Build the direct-space metric tensor and its inverse.
% The inverse metric is used to compute d spacings from HKL.
G = cell_metric_local(cif.a, cif.b, cif.c, cif.alpha, cif.beta, cif.gamma);
Gstar = inv(G);

% Resolve the requested output grids and determine the largest Q that
% must be represented.
[twoThetaGrid_deg, qGrid_Ainv, qMaxWanted] = resolve_grids_local(lambda_A, opts);

% The minimum d spacing needed is set by the largest desired Q:
%   Q = 2*pi/d  ->  d = 2*pi/Q
dMin = 2*pi / qMaxWanted;

% Choose a conservative HKL search box.
[hMax, kMax, lMax] = choose_hkl_box_local(Gstar, dMin, opts);

% Enumerate all HKL values in the search box.
[H,K,L] = ndgrid(-hMax:hMax, -kMax:kMax, -lMax:lMax);
H = H(:);
K = K(:);
L = L(:);

% Remove the forbidden zero reflection.
keep = ~(H == 0 & K == 0 & L == 0);
H = H(keep);
K = K(keep);
L = L(keep);

% Use the reciprocal-space metric to compute d spacings.
% For each HKL:
%   1/d^2 = [h k l] * Gstar * [h k l]'
HKL = [H K L];
invd2 = sum((HKL * Gstar) .* HKL, 2);
d_A = 1 ./ sqrt(invd2);

% Convert d spacing to Q.
Q_Ainv = 2*pi ./ d_A;

% Bragg condition:
%   sin(theta) = lambda / (2d)
% Only reflections with sin(theta) <= 1 are physically accessible.
arg = lambda_A ./ (2*d_A);

valid = (arg <= 1 + 1e-12) & ...
    (Q_Ainv >= min(qGrid_Ainv) - 1e-12) & ...
    (Q_Ainv <= qMaxWanted + 1e-12);

H = H(valid);
K = K(valid);
L = L(valid);
d_A = d_A(valid);
Q_Ainv = Q_Ainv(valid);
arg = arg(valid);

% Guard against tiny numerical overshoot above 1.
arg(arg > 1) = 1;

% Convert Bragg angle theta to diffraction angle 2theta.
twoTheta_deg = 2 * asind(arg);

nRef = numel(Q_Ainv);

% ---------------------------------------------------------------------
% Prepare atomic information for the structure-factor calculation.
% ---------------------------------------------------------------------
xyz = [atoms.x(:), atoms.y(:), atoms.z(:)];

if opts.corrections.useOccupancy
    occ = atoms.occ(:);
else
    occ = ones(size(atoms.occ(:)));
end

if opts.corrections.useDebyeWaller
    Biso = atoms.Biso(:);
else
    Biso = zeros(size(atoms.Biso(:)));
end

atomType = atoms.type(:);
[typeList, ~, typeID] = unique(atomType);

% Cache the CM coefficient record for each unique atom type so we do
% not repeat the lookup for every reflection and every atom.
cmCache = cell(numel(typeList), 1);
neutralList = cell(numel(typeList), 1);
Zlist = zeros(numel(typeList), 1);

for it = 1:numel(typeList)
    [cmCache{it}, neutralList{it}, Zlist(it)] = get_cm_record_local(typeList{it});
end

% Allocate structure factors and intensities.
Iraw = zeros(nRef, 1);
Fraw = zeros(nRef, 1);

nearEdgeWarned = false;

% ---------------------------------------------------------------------
% Main reflection loop: compute F(hkl) and |F(hkl)|^2
% ---------------------------------------------------------------------
for ii = 1:nRef
    h = H(ii);
    k = K(ii);
    l = L(ii);

    % Reduced scattering variable:
    %   s = sin(theta)/lambda = Q/(4*pi)
    s = Q_Ainv(ii) / (4*pi);

    % Complex phase factor from atomic positions.
    phase = exp(2*pi*1i*(h*xyz(:,1) + k*xyz(:,2) + l*xyz(:,3)));

    % Isotropic Debye-Waller factor.
    DW = exp(-Biso * s.^2);

    % Build the atomic form factor for every atom.
    % Base form:
    %   f = f0(Q)
    %
    % Optional near-edge form:
    %   f = f0(Q) + fp(E) + i*fpp(E)
    fj = zeros(size(occ));

    for it = 1:numel(typeList)
        idx = (typeID == it);

        % Non-dispersive X-ray form factor from CM coefficients.
        f0 = eval_cm_f0_local(cmCache{it}, s);

        fp = 0;
        fpp = 0;

        % Handle anomalous corrections if requested.
        if is_anomalous_enabled_local(opts.corrections.useAnomalous)
            if ~isempty(opts.corrections.anomalousFcn)
                [fp, fpp] = opts.corrections.anomalousFcn(neutralList{it}, EkeV);
            elseif strcmpi(opts.corrections.useAnomalous, 'auto')
                if is_near_edge_local(neutralList{it}, EkeV, opts.corrections.nearEdgeWindow_keV)
                    if opts.verbose && ~nearEdgeWarned
                        warning('cif2powder_CM:NearEdgeNoAnomalous', ...
                            'Photon energy appears near an absorption edge for at least one element, but no anomalous function was supplied. Using fp = 0 and fpp = 0.');
                        nearEdgeWarned = true;
                    end
                end
            end
        end

        fj(idx) = f0 + fp + 1i*fpp;
    end

    % Structure factor:
    %   F(hkl) = sum_j occ_j * DW_j * f_j * exp(phase_j)
    Fhkl = sum(occ .* DW .* fj .* phase);

    Fraw(ii) = Fhkl;
    Iraw(ii) = abs(Fhkl).^2;
end

% ---------------------------------------------------------------------
% Optional diffraction corrections applied to intensity.
% ---------------------------------------------------------------------
corrFactor = ones(size(Iraw));
theta_deg = twoTheta_deg / 2;

if opts.corrections.useLorentz
    corrFactor = corrFactor .* ...
        (1 ./ max(sind(theta_deg).^2 .* cosd(theta_deg), eps));
end

if opts.corrections.usePolarization
    corrFactor = corrFactor .* ((1 + cosd(twoTheta_deg).^2) / 2);
end

Icorr = Iraw .* corrFactor;

% ---------------------------------------------------------------------
% Merge nearly degenerate powder lines by Q.
% ---------------------------------------------------------------------
qTol = max(opts.mergeQtol_Ainv, 1e-6);
qKey = round(Q_Ainv / qTol);

[uKey, ~, ic] = unique(qKey); 
nLines = numel(uKey);

lineQ = accumarray(ic, Q_Ainv, [nLines 1], @mean);
lineTTH = accumarray(ic, twoTheta_deg, [nLines 1], @mean);
lineD = accumarray(ic, d_A, [nLines 1], @mean);
lineI = accumarray(ic, Icorr, [nLines 1], @sum);
lineMult = accumarray(ic, 1, [nLines 1], @sum);

hklRep = cell(nLines,1);
fwhmTTH_deg = zeros(nLines,1);
fwhmQ_Ainv = zeros(nLines,1);

% For each merged line, choose one representative HKL label and compute
% the local width for profile generation.
for jj = 1:nLines
    idx = find(ic == jj);
    [~, imax] = max(Icorr(idx));
    kk = idx(imax);

    hklRep{jj} = sprintf('(%d %d %d)', H(kk), K(kk), L(kk));

    fwhmTTH_deg(jj) = local_fwhm_tth_deg_local(lineTTH(jj), lambda_A, opts.profile);
    fwhmQ_Ainv(jj) = tth_fwhm_to_q_fwhm_local(lineTTH(jj), fwhmTTH_deg(jj), lambda_A);
end

% Sort merged lines by Q.
[lineQ, ord] = sort(lineQ);
lineTTH = lineTTH(ord);
lineD = lineD(ord);
lineI = lineI(ord);
lineMult = lineMult(ord);
hklRep = hklRep(ord);
fwhmTTH_deg = fwhmTTH_deg(ord);
fwhmQ_Ainv = fwhmQ_Ainv(ord);

% ---------------------------------------------------------------------
% Build the final broadened profiles.
% ---------------------------------------------------------------------
tthProfile = build_profile_variable_width_local( ...
    twoThetaGrid_deg(:), lineTTH, lineI, fwhmTTH_deg, opts.profile);

qProfile = build_profile_variable_width_local( ...
    qGrid_Ainv(:), lineQ, lineI, fwhmQ_Ainv, opts.profile);

% ---------------------------------------------------------------------
% Pack output.
% ---------------------------------------------------------------------
out = struct;
out.lambda_A = lambda_A;
out.energy_keV = EkeV;
out.options = opts;
out.cell = rmfield(cif, {'symops','atom'});
out.atomsExpanded = atoms;

out.linesRaw = struct;
out.linesRaw.h = H;
out.linesRaw.k = K;
out.linesRaw.l = L;
out.linesRaw.d_A = d_A;
out.linesRaw.Q_Ainv = Q_Ainv;
out.linesRaw.twoTheta_deg = twoTheta_deg;
out.linesRaw.F = Fraw;
out.linesRaw.I = Icorr;

out.linesMerged = struct;
out.linesMerged.hkl_rep = hklRep;
out.linesMerged.d_A = lineD;
out.linesMerged.Q_Ainv = lineQ;
out.linesMerged.twoTheta_deg = lineTTH;
out.linesMerged.I = lineI;
out.linesMerged.multiplicity = lineMult;
out.linesMerged.fwhmTwoTheta_deg = fwhmTTH_deg;
out.linesMerged.fwhmQ_Ainv = fwhmQ_Ainv;

out.tthGrid_deg = twoThetaGrid_deg(:);
out.tthProfile = tthProfile(:);
out.tthProfileNorm = tthProfile(:) / max(tthProfile(:) + eps);

out.qGrid_Ainv = qGrid_Ainv(:);
out.qProfile = qProfile(:);
out.qProfileNorm = qProfile(:) / max(qProfile(:) + eps);
out.qPeaks = out.linesMerged.Q_Ainv(:);
out.qGrid=qGrid_Ainv(:);


end


% =========================================================================
% DEFAULT OPTIONS
% =========================================================================
function opts = set_default_opts_local(opts)

if ~isfield(opts, 'twoThetaGrid_deg') || isempty(opts.twoThetaGrid_deg)
    opts.twoThetaGrid_deg = 0.1:0.02:90;
end

if ~isfield(opts, 'qGrid_Ainv') || isempty(opts.qGrid_Ainv)
    opts.qGrid_Ainv = [];
end

if ~isfield(opts, 'mergeQtol_Ainv') || isempty(opts.mergeQtol_Ainv)
    opts.mergeQtol_Ainv = 1e-4;
end

if ~isfield(opts, 'hklMax') || isempty(opts.hklMax)
    opts.hklMax = [];
end

if ~isfield(opts, 'verbose') || isempty(opts.verbose)
    opts.verbose = true;
end

if ~isfield(opts, 'corrections') || isempty(opts.corrections)
    opts.corrections = struct;
end

if ~isfield(opts.corrections, 'useOccupancy') || isempty(opts.corrections.useOccupancy)
    opts.corrections.useOccupancy = true;
end

if ~isfield(opts.corrections, 'useDebyeWaller') || isempty(opts.corrections.useDebyeWaller)
    opts.corrections.useDebyeWaller = true;
end

if ~isfield(opts.corrections, 'useLorentz') || isempty(opts.corrections.useLorentz)
    opts.corrections.useLorentz = true;
end

if ~isfield(opts.corrections, 'usePolarization') || isempty(opts.corrections.usePolarization)
    opts.corrections.usePolarization = false;
end

if ~isfield(opts.corrections, 'useAnomalous') || isempty(opts.corrections.useAnomalous)
    opts.corrections.useAnomalous = 'auto';
end

if ~isfield(opts.corrections, 'anomalousFcn') || isempty(opts.corrections.anomalousFcn)
    opts.corrections.anomalousFcn = [];
end

if ~isfield(opts.corrections, 'nearEdgeWindow_keV') || isempty(opts.corrections.nearEdgeWindow_keV)
    opts.corrections.nearEdgeWindow_keV = 1.0;
end

if ~isfield(opts, 'profile') || isempty(opts.profile)
    opts.profile = struct;
end

if ~isfield(opts.profile, 'widthModel') || isempty(opts.profile.widthModel)
    opts.profile.widthModel = 'fixed';
end

if ~isfield(opts.profile, 'shape') || isempty(opts.profile.shape)
    opts.profile.shape = 'gaussian';
end

if ~isfield(opts.profile, 'fixedFWHM_deg') || isempty(opts.profile.fixedFWHM_deg)
    opts.profile.fixedFWHM_deg = 0.06;
end

if ~isfield(opts.profile, 'U') || isempty(opts.profile.U)
    opts.profile.U = 0.002;
end

if ~isfield(opts.profile, 'V') || isempty(opts.profile.V)
    opts.profile.V = 0;
end

if ~isfield(opts.profile, 'W') || isempty(opts.profile.W)
    opts.profile.W = 0.003;
end

if ~isfield(opts.profile, 'useSize') || isempty(opts.profile.useSize)
    opts.profile.useSize = false;
end

if ~isfield(opts.profile, 'size_A') || isempty(opts.profile.size_A)
    opts.profile.size_A = Inf;
end

if ~isfield(opts.profile, 'shapeK') || isempty(opts.profile.shapeK)
    opts.profile.shapeK = 0.9;
end

if ~isfield(opts.profile, 'useStrain') || isempty(opts.profile.useStrain)
    opts.profile.useStrain = false;
end

if ~isfield(opts.profile, 'microstrain') || isempty(opts.profile.microstrain)
    opts.profile.microstrain = 0;
end

if ~isfield(opts.profile, 'eta') || isempty(opts.profile.eta)
    opts.profile.eta = 0.3;
end

opts.profile.eta = max(0, min(1, opts.profile.eta));
end


% =========================================================================
% GRID RESOLUTION
% =========================================================================
function [twoThetaGrid_deg, qGrid_Ainv, qMaxWanted] = resolve_grids_local(lambda_A, opts)

twoThetaGrid_deg = opts.twoThetaGrid_deg(:).';
qGrid_Ainv = opts.qGrid_Ainv(:).';

qMaxByEnergy = 4*pi / lambda_A;

if isempty(twoThetaGrid_deg) && isempty(qGrid_Ainv)
    twoThetaGrid_deg = 0.1:0.02:90;
end

% If no Q grid is given, derive one that spans the requested 2theta range.
if isempty(qGrid_Ainv)
    qFromTTH = 4*pi*sind(twoThetaGrid_deg/2) / lambda_A;
    qGrid_Ainv = linspace(min(qFromTTH), min(max(qFromTTH), qMaxByEnergy), numel(twoThetaGrid_deg));
end

% If no 2theta grid is given, derive one from the requested Q grid.
if isempty(twoThetaGrid_deg)
    arg = qGrid_Ainv * lambda_A / (4*pi);
    arg(arg > 1) = NaN;
    twoThetaGrid_deg = 2*asind(arg);
    twoThetaGrid_deg = twoThetaGrid_deg(isfinite(twoThetaGrid_deg));
end

qMaxFromTTH = 4*pi*sind(max(twoThetaGrid_deg)/2) / lambda_A;
qMaxWanted = min(max([qGrid_Ainv(:); qMaxFromTTH]), qMaxByEnergy);

qGrid_Ainv = qGrid_Ainv(qGrid_Ainv <= qMaxByEnergy + 1e-12);
end


% =========================================================================
% HKL SEARCH BOX
% =========================================================================
function [hMax, kMax, lMax] = choose_hkl_box_local(Gstar, dMin, opts)

if ~isempty(opts.hklMax)
    if isscalar(opts.hklMax)
        hMax = opts.hklMax;
        kMax = opts.hklMax;
        lMax = opts.hklMax;
    else
        hMax = opts.hklMax(1);
        kMax = opts.hklMax(2);
        lMax = opts.hklMax(3);
    end
    return
end

% Reciprocal basis lengths from diagonal elements of Gstar.
astar = sqrt(Gstar(1,1));
bstar = sqrt(Gstar(2,2));
cstar = sqrt(Gstar(3,3));

% Conservative rectangular HKL search box.
% The exact admissible reflections are selected afterward using:
%
%   [h k l] * Gstar * [h k l]' <= 1/dMin^2
%
% The rectangular box is intentionally loose but safe.
hMax = max(1, ceil(1/(dMin*astar)) + 1);
kMax = max(1, ceil(1/(dMin*bstar)) + 1);
lMax = max(1, ceil(1/(dMin*cstar)) + 1);
end


% =========================================================================
% LINE-WIDTH MODEL
% =========================================================================
function fwhm_deg = local_fwhm_tth_deg_local(twoTheta_deg, lambda_A, profile)

theta_rad = (twoTheta_deg/2) * pi/180;

switch lower(profile.widthModel)
    case 'fixed'
        % Every peak has the same width in 2theta.
        fwhm_deg = profile.fixedFWHM_deg;

    case 'caglioti'
        % Standard empirical instrumental-width model:
        %   FWHM^2 = U*tan(theta)^2 + V*tan(theta) + W
        tanth = tan(theta_rad);
        H2 = profile.U * tanth.^2 + profile.V * tanth + profile.W;
        H2 = max(H2, eps);
        fwhm_deg = sqrt(H2);

    otherwise
        error('cif2powder_CM:BadWidthModel', ...
            'Unknown profile.widthModel: %s', profile.widthModel);
end

% Optional Scherrer-type size broadening.
if profile.useSize && isfinite(profile.size_A) && profile.size_A > 0
    betaSize_rad = profile.shapeK * lambda_A ./ max(profile.size_A * cos(theta_rad), eps);
    betaSize_deg = betaSize_rad * 180/pi;
    fwhm_deg = sqrt(fwhm_deg.^2 + betaSize_deg.^2);
end

% Optional microstrain broadening.
if profile.useStrain && profile.microstrain > 0
    betaStrain_rad = 4 * profile.microstrain * tan(theta_rad);
    betaStrain_deg = betaStrain_rad * 180/pi;
    fwhm_deg = sqrt(fwhm_deg.^2 + betaStrain_deg.^2);
end

fwhm_deg = max(fwhm_deg, 1e-6);
end

function fwhmQ = tth_fwhm_to_q_fwhm_local(twoTheta_deg, fwhmTTH_deg, lambda_A)
% Convert a local width in 2theta to a local width in Q using the
% derivative dQ/d(2theta).
theta_deg = twoTheta_deg / 2;
dQd2theta_deg = (2*pi/lambda_A) * cosd(theta_deg) * (pi/180);
fwhmQ = abs(dQd2theta_deg) * fwhmTTH_deg;
fwhmQ = max(fwhmQ, 1e-8);
end


% =========================================================================
% PROFILE BUILDER
% =========================================================================
function prof = build_profile_variable_width_local(xGrid, x0, amp, fwhm, profile)

prof = zeros(size(xGrid));

for ii = 1:numel(x0)
    switch lower(profile.shape)
        case 'gaussian'
            prof = prof + amp(ii) * gaussian_unit_area_local(xGrid, x0(ii), fwhm(ii));

        case 'pvoigt'
            eta = profile.eta;
            prof = prof + amp(ii) * ...
                ((1-eta) * gaussian_unit_area_local(xGrid, x0(ii), fwhm(ii)) + ...
                eta    * lorentz_unit_area_local(xGrid, x0(ii), fwhm(ii)));

        otherwise
            error('cif2powder_CM:BadProfileShape', ...
                'Unknown profile.shape: %s', profile.shape);
    end
end
end

function y = gaussian_unit_area_local(x, x0, fwhm)
sigma = fwhm / (2*sqrt(2*log(2)));
y = exp(-0.5*((x-x0)/sigma).^2) / (sigma * sqrt(2*pi));
end

function y = lorentz_unit_area_local(x, x0, fwhm)
gamma = fwhm / 2;
y = (gamma/pi) ./ ((x-x0).^2 + gamma.^2);
end


% =========================================================================
% CM COEFFICIENT ACCESS
% =========================================================================
function [rec, neutralName, Z] = get_cm_record_local(typeName)

if isstring(typeName)
    typeName = char(typeName);
end
typeName = strtrim(typeName);

exactName = cleanup_species_name_local(typeName);

% Extract a neutral element symbol from a CIF atom type, for example in
% the LaB6 case, B0 -> B, La0 -> La, La3+ -> La
tok = regexp(exactName, '^[A-Z][a-z]?', 'match', 'once');
if isempty(tok)
    error('cif2powder_CM:BadSpecies', ...
        'Could not infer element symbol from "%s".', typeName);
end
neutralName = tok;

% Try the exact name first.
try
    [rec, Z] = CMcoef_cif(exactName);
    validate_cm_output_local(rec, Z, exactName);
    rec = rec(:);
    return
catch
end

% Fall back to the neutral element symbol.
try
    [rec, Z] = CMcoef_cif(neutralName);
    validate_cm_output_local(rec, Z, neutralName);
    rec = rec(:);
    return
catch ME
    error('cif2powder_CM:NoCMMatch', ...
        'Could not find CM coefficients for CIF type "%s" (cleaned "%s", neutral "%s"). Original error: %s', ...
        typeName, exactName, neutralName, ME.message);
end
end

function validate_cm_output_local(rec, Z, nameUsed)
if isempty(rec) || numel(rec) ~= 11
    error('cif2powder_CM:BadCMcoef', ...
        'CMcoef_cif returned an invalid coefficient vector for "%s".', nameUsed);
end
if isempty(Z) || ~isscalar(Z) || ~isfinite(Z)
    error('cif2powder_CM:BadZ', ...
        'CMcoef_cif returned an invalid Z for "%s".', nameUsed);
end
end

function f0 = eval_cm_f0_local(rec, s)
% CM record layout:
%   [a1 a2 a3 a4 a5 c b1 b2 b3 b4 b5]'
a = rec(1:5);
c = rec(6);
b = rec(7:11);

f0 = c + sum(a .* exp(-b * s.^2));
end

function tf = is_anomalous_enabled_local(flag)
if islogical(flag)
    tf = flag;
elseif ischar(flag)
    tf = strcmpi(flag, 'auto') || strcmpi(flag, 'true') || strcmpi(flag, 'yes');
else
    tf = false;
end
end


% =========================================================================
% APPROXIMATE EDGE CHECKER
% =========================================================================
function tf = is_near_edge_local(elemNeutral, EkeV, window_keV)

tf = false;
edges = approximate_edge_energies_local(elemNeutral);

if isempty(edges)
    return
end

if any(abs(edges(:) - EkeV) <= window_keV)
    tf = true;
end
end

function edges_keV = approximate_edge_energies_local(elemNeutral)
% Minimal local list, only for warning logic in "auto" mode.
% These values are NOT a substitute for a true fp/fpp table.
switch elemNeutral
    case 'B'
        edges_keV = 0.188;
    case 'La'
        edges_keV = [5.48 5.89 6.27 38.9];
    otherwise
        edges_keV = [];
end
end


% =========================================================================
% CIF PARSING
% =========================================================================
function cif = read_cif_basic_local(cifFile)

txt = fileread(cifFile);
txt = strrep(txt, sprintf('\r\n'), newline);
txt = strrep(txt, sprintf('\r'), newline);

lines = regexp(txt, '\n', 'split');
lines = lines(:);

% Remove trailing comments beginning with #.
for i = 1:numel(lines)
    k = strfind(lines{i}, '#');
    if ~isempty(k)
        lines{i} = lines{i}(1:k(1)-1);
    end
end

cif = struct;
cif.a = get_tag_number_local(lines, {'_cell_length_a'});
cif.b = get_tag_number_local(lines, {'_cell_length_b'});
cif.c = get_tag_number_local(lines, {'_cell_length_c'});
cif.alpha = get_tag_number_local(lines, {'_cell_angle_alpha'});
cif.beta = get_tag_number_local(lines, {'_cell_angle_beta'});
cif.gamma = get_tag_number_local(lines, {'_cell_angle_gamma'});

cif.symops = {};
cif.atom = struct('label', {{}}, 'type', {{}}, 'x', [], 'y', [], 'z', [], 'occ', [], 'Biso', []);

i = 1;
while i <= numel(lines)
    s = strtrim(lines{i});

    if strcmpi(s, 'loop_')
        j = i + 1;
        tags = {};

        % Read the loop column names.
        while j <= numel(lines)
            sj = strtrim(lines{j});
            if ~isempty(sj) && sj(1) == '_'
                tok = tokenize_cif_line_local(sj);
                tags{end+1} = tok{1}; %#ok<AGROW>
                j = j + 1;
            else
                break
            end
        end

        % Read the loop rows until the next loop or tag block.
        rows = {};
        while j <= numel(lines)
            sj = strtrim(lines{j});
            if isempty(sj)
                j = j + 1;
                continue
            end

            if strcmpi(sj, 'loop_') || ...
                    (~isempty(sj) && sj(1) == '_') || ...
                    strncmpi(sj, 'data_', 5)
                break
            end

            rows{end+1} = tokenize_cif_line_local(lines{j}); %#ok<AGROW>
            j = j + 1;
        end

        % Look for symmetry operation loops.
        iSym = find(strcmpi(tags, '_space_group_symop_operation_xyz') | ...
            strcmpi(tags, '_space_group_symop.operation_xyz') | ...
            strcmpi(tags, '_symmetry_equiv_pos_as_xyz'), 1);

        if ~isempty(iSym)
            symops = {};
            for r = 1:numel(rows)
                row = rows{r};
                if numel(row) >= iSym
                    symops{end+1} = strip_quotes_local(row{iSym}); %#ok<AGROW>
                end
            end
            cif.symops = symops;
        end

        % Look for atomic position loops.
        hasX = any(strcmpi(tags, '_atom_site_fract_x'));
        hasY = any(strcmpi(tags, '_atom_site_fract_y'));
        hasZ = any(strcmpi(tags, '_atom_site_fract_z'));

        if hasX && hasY && hasZ
            ix = find(strcmpi(tags, '_atom_site_fract_x'), 1);
            iy = find(strcmpi(tags, '_atom_site_fract_y'), 1);
            iz = find(strcmpi(tags, '_atom_site_fract_z'), 1);

            ilab = find(strcmpi(tags, '_atom_site_label'), 1);
            ityp = find(strcmpi(tags, '_atom_site_type_symbol'), 1);
            iocc = find(strcmpi(tags, '_atom_site_occupancy'), 1);

            iB = find(strcmpi(tags, '_atom_site_B_iso_or_equiv'), 1);
            iU = find(strcmpi(tags, '_atom_site_U_iso_or_equiv'), 1);

            label = {};
            type = {};
            x = [];
            y = [];
            z = [];
            occ = [];
            Biso = [];

            for r = 1:numel(rows)
                row = rows{r};

                if numel(row) < max([ix iy iz])
                    continue
                end

                xr = parse_cif_number_local(row{ix});
                yr = parse_cif_number_local(row{iy});
                zr = parse_cif_number_local(row{iz});

                if any(isnan([xr yr zr]))
                    continue
                end

                if ~isempty(ilab) && numel(row) >= ilab
                    lab = strip_quotes_local(row{ilab});
                else
                    lab = sprintf('A%d', r);
                end

                if ~isempty(ityp) && numel(row) >= ityp
                    typ = strip_quotes_local(row{ityp});
                else
                    typ = infer_symbol_from_label_local(lab);
                end

                if ~isempty(iocc) && numel(row) >= iocc
                    occr = parse_cif_number_local(row{iocc});
                    if isnan(occr)
                        occr = 1;
                    end
                else
                    occr = 1;
                end

                if ~isempty(iB) && numel(row) >= iB
                    Br = parse_cif_number_local(row{iB});
                    if isnan(Br)
                        Br = 0;
                    end
                elseif ~isempty(iU) && numel(row) >= iU
                    Ur = parse_cif_number_local(row{iU});
                    if isnan(Ur)
                        Ur = 0;
                    end
                    Br = 8*pi^2*Ur;
                else
                    Br = 0;
                end

                label{end+1,1} = lab; %#ok<AGROW>
                type{end+1,1} = typ; %#ok<AGROW>
                x(end+1,1) = xr; %#ok<AGROW>
                y(end+1,1) = yr; %#ok<AGROW>
                z(end+1,1) = zr; %#ok<AGROW>
                occ(end+1,1) = occr; %#ok<AGROW>
                Biso(end+1,1) = Br; %#ok<AGROW>
            end

            cif.atom.label = label;
            cif.atom.type = type;
            cif.atom.x = x;
            cif.atom.y = y;
            cif.atom.z = z;
            cif.atom.occ = occ;
            cif.atom.Biso = Biso;
        end

        i = j;
    else
        i = i + 1;
    end
end

if isempty(cif.symops)
    cif.symops = {'x,y,z'};
end

if isempty(cif.atom.x)
    error('cif2powder_CM:NoAtoms', ...
        'No atomic fractional coordinates were found in the CIF.');
end
end

function toks = tokenize_cif_line_local(s)
toks = regexp(s, '''[^'']*''|"[^"]*"|\S+', 'match');
end

function s = strip_quotes_local(s)
s = strtrim(s);
if isempty(s)
    return
end

if (s(1) == '''' && s(end) == '''') || (s(1) == '"' && s(end) == '"')
    s = s(2:end-1);
end
end

function x = parse_cif_number_local(s)
s = strip_quotes_local(strtrim(s));

if isempty(s) || strcmp(s,'.') || strcmp(s,'?')
    x = NaN;
    return
end

% Remove crystallographic uncertainty notation, for example 4.1(2)->4.1
s = regexprep(s, '\([^\)]*\)', '');
x = str2double(s);
end

function val = get_tag_number_local(lines, tagList)
val = NaN;

for it = 1:numel(tagList)
    tag = tagList{it};

    for i = 1:numel(lines)
        s = strtrim(lines{i});

        if numel(s) >= numel(tag) && strcmpi(s(1:numel(tag)), tag)
            tok = tokenize_cif_line_local(s);

            if numel(tok) >= 2
                val = parse_cif_number_local(tok{2});
            elseif i < numel(lines)
                val = parse_cif_number_local(strtrim(lines{i+1}));
            end

            if ~isnan(val)
                return
            end
        end
    end
end

if isnan(val)
    error('cif2powder_CM:MissingTag', ...
        'Missing CIF tag: %s', tagList{1});
end
end


% =========================================================================
% ASYMMETRIC-UNIT EXPANSION
% =========================================================================
function atomsFull = expand_asymmetric_unit_local(cif)

tol = 1e-6;

xAll = [];
yAll = [];
zAll = [];
occAll = [];
BisoAll = [];
typeAll = {};
labelAll = {};

for ia = 1:numel(cif.atom.x)
    p0 = [cif.atom.x(ia), cif.atom.y(ia), cif.atom.z(ia)];
    xyz = zeros(numel(cif.symops), 3);

    for is = 1:numel(cif.symops)
        xyz(is,:) = apply_symop_local(cif.symops{is}, p0);
    end

    % Wrap all generated positions back into the unit cell.
    xyz = mod(xyz, 1);

    % Remove duplicates generated by symmetry or special positions.
    xyzKey = round(xyz / tol) * tol;
    [~, iu] = unique(xyzKey, 'rows', 'stable');
    xyz = xyz(iu,:);

    nEq = size(xyz,1);

    xAll = [xAll; xyz(:,1)];
    yAll = [yAll; xyz(:,2)];
    zAll = [zAll; xyz(:,3)];
    occAll = [occAll; cif.atom.occ(ia) * ones(nEq,1)];
    BisoAll = [BisoAll; cif.atom.Biso(ia) * ones(nEq,1)];
    typeAll = [typeAll; repmat(cif.atom.type(ia), nEq, 1)];
    labelAll = [labelAll; repmat(cif.atom.label(ia), nEq, 1)];
end

atomsFull = struct;
atomsFull.x = xAll;
atomsFull.y = yAll;
atomsFull.z = zAll;
atomsFull.occ = occAll;
atomsFull.Biso = BisoAll;
atomsFull.type = typeAll;
atomsFull.label = labelAll;
end

function xyzOut = apply_symop_local(symStr, xyzIn)

s = lower(strtrim(symStr));
s = strrep(s, '''', '');
s = strrep(s, '"', '');
s = strrep(s, ' ', '');

parts = strsplit(s, ',');
if numel(parts) ~= 3
    error('cif2powder_CM:BadSymOp', ...
        'Could not parse symmetry operation: %s', symStr);
end

x = xyzIn(1);
y = xyzIn(2);
z = xyzIn(3);

xyzOut = zeros(1,3);
xyzOut(1) = eval_sym_expr_local(parts{1}, x, y, z);
xyzOut(2) = eval_sym_expr_local(parts{2}, x, y, z);
xyzOut(3) = eval_sym_expr_local(parts{3}, x, y, z);
end

function val = eval_sym_expr_local(expr, x, y, z)

expr = strrep(expr, '-', '+-');
if ~isempty(expr) && expr(1) == '+'
    expr = expr(2:end);
end

toks = strsplit(expr, '+');
val = 0;

for it = 1:numel(toks)
    tok = toks{it};

    if isempty(tok)
        continue
    end

    if contains(tok, 'x')
        c = strrep(tok, 'x', '');
        coeff = parse_linear_coeff_local(c);
        val = val + coeff * x;

    elseif contains(tok, 'y')
        c = strrep(tok, 'y', '');
        coeff = parse_linear_coeff_local(c);
        val = val + coeff * y;

    elseif contains(tok, 'z')
        c = strrep(tok, 'z', '');
        coeff = parse_linear_coeff_local(c);
        val = val + coeff * z;

    else
        val = val + parse_fraction_local(tok);
    end
end
end

function coeff = parse_linear_coeff_local(s)
if isempty(s)
    coeff = 1;
elseif strcmp(s, '-')
    coeff = -1;
else
    coeff = str2double(s);
    if isnan(coeff)
        if contains(s, '-')
            coeff = -1;
        else
            coeff = 1;
        end
    end
end
end

function num = parse_fraction_local(s)
s = strtrim(s);

if isempty(s)
    num = 0;
    return
end

if contains(s, '/')
    p = strsplit(s, '/');
    num = str2double(p{1}) / str2double(p{2});
else
    num = str2double(s);
end

if isnan(num)
    error('cif2powder_CM:BadNumericToken', ...
        'Could not parse numeric token: %s', s);
end
end


% =========================================================================
% METRIC TENSOR
% =========================================================================
function G = cell_metric_local(a,b,c,alpha,beta,gamma)

ca = cosd(alpha);
cb = cosd(beta);
cg = cosd(gamma);

G = [a^2,    a*b*cg, a*c*cb; ...
    a*b*cg, b^2,    b*c*ca; ...
    a*c*cb, b*c*ca, c^2];
end


% =========================================================================
% SYMBOL CLEANUP
% =========================================================================
function s = cleanup_species_name_local(s)
if isstring(s)
    s = char(s);
end
s = strtrim(s);
s = regexprep(s, '\s+', '');
s = strrep(s, '"', '');
s = strrep(s, '''', '');
end

function sym = infer_symbol_from_label_local(lab)
tok = regexp(lab, '[A-Z][a-z]?', 'match', 'once');
if isempty(tok)
    error('cif2powder_CM:BadLabelSymbol', ...
        'Could not infer chemical symbol from label "%s".', lab);
end
sym = tok;
end


function [f, varargout] = CMcoef_cif(name)
% Retrieve the Cromer-Mann-Waasmaier-Kirfel coefficients
% The data is from f0_WaasKirf.dat ( ftp.esrf.eu/pub/scisoft/DabaxFiles/ )
% and was computed for neutral atoms and ions for 0.0 to 6.0 A-1.
% The non-dispersive part of the scattering factor is approximated by:
%
%    f0[k] = c + [SUM a_i*EXP(-b_i*(k^2)) ]
%                i=1,5
%
% where k = sin(theta)/lambda = q/4pi, and c, a_i and b_i
% are the coefficients tabulated in T in the order:
% [AtomName Z a1  a2  a3  a4  a5  c  b1  b2  b3  b4  b5]
%
% Inputs:
%  name   - the name of an atom  (char)
%
% Outputs:
%
%  f      - the coef vector (a1..a5 c b1..b5)
%  z      - the Z number of the atom (optional)
%  Adi Natan


T ={'H'	1	0.413048000000000	0.294953000000000	0.187491000000000	0.0807010000000000	0.0237360000000000	4.90000000000000e-05	15.5699460000000	32.3984680000000	5.71140400000000	61.8898740000000	1.33411800000000
    'H1-'	1	0.702260000000000	0.763666000000000	0.248678000000000	0.261323000000000	0.0230170000000000	0.000425000000000000	23.9456040000000	74.8979190000000	6.77328900000000	233.583450000000	1.33753100000000
    'He'	2	0.732354000000000	0.753896000000000	0.283819000000000	0.190003000000000	0.0391390000000000	0.000487000000000000	11.5539180000000	4.59583100000000	1.54629900000000	26.4639640000000	0.377523000000000
    'Li'	3	0.974637000000000	0.158472000000000	0.811855000000000	0.262416000000000	0.790108000000000	0.00254200000000000	4.33494600000000	0.342451000000000	97.1029660000000	201.363831000000	1.40923400000000
    'Li1+'	3	0.432724000000000	0.549257000000000	0.376575000000000	-0.336481000000000	0.976060000000000	0.00176400000000000	0.260367000000000	1.04283600000000	7.88529400000000	0.260368000000000	3.04253900000000
    'Be'	4	1.53371200000000	0.638283000000000	0.601052000000000	0.106139000000000	1.11841400000000	0.00251100000000000	42.6620790000000	0.595420000000000	99.1064990000000	0.151340000000000	1.84309300000000
    'Be2+'	4	3.05543000000000	-2.37261700000000	1.04491400000000	0.544233000000000	0.381737000000000	-0.653773000000000	0.00122600000000000	0.00122700000000000	1.54210600000000	0.456279000000000	4.04747900000000
    'B'	    5	2.08518500000000	1.06458000000000	1.06278800000000	0.140515000000000	0.641784000000000	0.00382300000000000	23.4940680000000	1.13789400000000	61.2389760000000	0.114886000000000	0.399036000000000
    'C'	    6	2.65750600000000	1.07807900000000	1.49090900000000	-4.24107000000000	0.713791000000000	4.29798300000000	14.7807580000000	0.776775000000000	42.0868420000000	-0.000294000000000000	0.239535000000000
    'Cval'	6	1.25848900000000	0.728215000000000	1.11985600000000	2.16813300000000	0.705239000000000	0.0197220000000000	10.6837690000000	0.208177000000000	0.836097000000000	24.6037040000000	58.9542730000000
    'N'	    7	11.8937800000000	3.27747900000000	1.85809200000000	0.858927000000000	0.912985000000000	-11.8049020000000	0.000158000000000000	10.2327230000000	30.3446900000000	0.656065000000000	0.217287000000000
    'O'	    8	2.96042700000000	2.50881800000000	0.637853000000000	0.722838000000000	1.14275600000000	0.0270140000000000	14.1822590000000	5.93685800000000	0.112726000000000	34.9584810000000	0.390240000000000
    'O1-'	8	3.10693400000000	3.23514200000000	1.14888600000000	0.783981000000000	0.676953000000000	0.0461360000000000	19.8680800000000	6.96025200000000	0.170043000000000	65.6935120000000	0.630757000000000
    'O2-'	8	3.99024700000000	2.30056300000000	0.607200000000000	1.90788200000000	1.16708000000000	0.0254290000000000	16.6399560000000	5.63681900000000	0.108493000000000	47.2997090000000	0.379984000000000
    'F'	    9	3.51194300000000	2.77224400000000	0.678385000000000	0.915159000000000	1.08926100000000	0.0325570000000000	10.6878590000000	4.38046600000000	0.0939820000000000	27.2552030000000	0.313066000000000
    'F1-'	9	0.457649000000000	3.84156100000000	1.43277100000000	0.801876000000000	3.39504100000000	0.0695250000000000	0.917243000000000	5.50780300000000	0.164955000000000	51.0762060000000	15.8216790000000
    'Ne'	10	4.18374900000000	2.90572600000000	0.520513000000000	1.13564100000000	1.22806500000000	0.0255760000000000	8.17545700000000	3.25253600000000	0.0632950000000000	21.8139100000000	0.224952000000000
    'Na'	11	4.91012700000000	3.08178300000000	1.26206700000000	1.09893800000000	0.560991000000000	0.0797120000000000	3.28143400000000	9.11917800000000	0.102763000000000	132.013947000000	0.405878000000000
    'Na1+'	11	3.14869000000000	4.07398900000000	0.767888000000000	0.995612000000000	0.968249000000000	0.0453000000000000	2.59498700000000	6.04692500000000	0.0701390000000000	14.1226570000000	0.217037000000000
    'Mg'	12	4.70897100000000	1.19481400000000	1.55815700000000	1.17041300000000	3.23940300000000	0.126842000000000	4.87520700000000	108.506081000000	0.111516000000000	48.2924080000000	1.92817100000000
    'Mg2+'	12	3.06291800000000	4.13510600000000	0.853742000000000	1.03679200000000	0.852520000000000	0.0588510000000000	2.01580300000000	4.41794100000000	0.0653070000000000	9.66971000000000	0.187818000000000
    'Al'	13	4.73079600000000	2.31395100000000	1.54198000000000	1.11756400000000	3.15475400000000	0.139509000000000	3.62893100000000	43.0511670000000	0.0959600000000000	108.932388000000	1.55591800000000
    'Al3+'	13	4.13201500000000	0.912049000000000	1.10242500000000	0.614876000000000	3.21913600000000	0.0193970000000000	3.52864100000000	7.37834400000000	0.133708000000000	0.0390650000000000	1.64472800000000
    'Si'	14	5.27532900000000	3.19103800000000	1.51151400000000	1.35684900000000	2.51911400000000	0.145073000000000	2.63133800000000	33.7307280000000	0.0811190000000000	86.2886430000000	1.17008700000000
    'Siva'	14	2.87903300000000	3.07296000000000	1.51598100000000	1.39003000000000	4.99505100000000	0.146030000000000	1.23971300000000	38.7062760000000	0.0814810000000000	93.6163330000000	2.77029300000000
    'Si4+'	14	3.67672200000000	3.82849600000000	1.25803300000000	0.419024000000000	0.720421000000000	0.0972660000000000	1.44685100000000	3.01314400000000	0.0643970000000000	0.206254000000000	5.97022200000000
    'P'	    15	1.95054100000000	4.14693000000000	1.49456000000000	1.52204200000000	5.72971100000000	0.155233000000000	0.908139000000000	27.0449520000000	0.0712800000000000	67.5201870000000	1.98117300000000
    'S'	    16	6.37215700000000	5.15456800000000	1.47373200000000	1.63507300000000	1.20937200000000	0.154722000000000	1.51434700000000	22.0925270000000	0.0613730000000000	55.4451750000000	0.646925000000000
    'Cl'	17	1.44607100000000	6.87060900000000	6.15180100000000	1.75034700000000	0.634168000000000	0.146773000000000	0.0523570000000000	1.19316500000000	18.3434160000000	46.3983960000000	0.401005000000000
    'Cl1-'	17	1.06180200000000	7.13988600000000	6.52427100000000	2.35562600000000	35.8294030000000	-34.9166030000000	0.144727000000000	1.17179500000000	19.4676550000000	60.3203010000000	0.000436000000000000
    'Ar'	18	7.18800400000000	6.63845400000000	0.454180000000000	1.92959300000000	1.52365400000000	0.265954000000000	0.956221000000000	15.3398770000000	15.3398620000000	39.0438230000000	0.0624090000000000
    'K'	    19	8.16399100000000	7.14694500000000	1.07014000000000	0.877316000000000	1.48643400000000	0.253614000000000	12.8163230000000	0.808945000000000	210.327011000000	39.5976520000000	0.0528210000000000
    'K1+'	19	-17.6093390000000	1.49487300000000	7.15030500000000	10.8995690000000	15.8082280000000	0.257164000000000	18.8409790000000	0.0534530000000000	0.812940000000000	22.2641050000000	14.3515930000000
    'Ca'	20	8.59365500000000	1.47732400000000	1.43625400000000	1.18283900000000	7.11325800000000	0.196255000000000	10.4606440000000	0.0418910000000000	81.3903810000000	169.847839000000	0.688098000000000
    'Ca2+'	20	8.50144100000000	12.8804830000000	9.76509500000000	7.15666900000000	0.711160000000000	-21.0131870000000	10.5258480000000	-0.00403300000000000	0.0106920000000000	0.684443000000000	27.2317710000000
    'Sc'	21	1.47656600000000	1.48727800000000	1.60018700000000	9.17746300000000	7.09975000000000	0.157765000000000	53.1310230000000	0.0353250000000000	137.319489000000	9.09803100000000	0.602102000000000
    'Sc3+'	21	7.10434800000000	1.51148800000000	-53.6697730000000	38.4048160000000	24.5322400000000	0.118642000000000	0.601957000000000	0.0333860000000000	12.5721380000000	10.8597360000000	14.1252300000000
    'Ti'	22	9.81852400000000	1.52264600000000	1.70310100000000	1.76877400000000	7.08255500000000	0.102473000000000	8.00187900000000	0.0297630000000000	39.8854220000000	120.157997000000	0.532405000000000
    'Ti2+'	22	7.04011900000000	1.49628500000000	9.65730400000000	0.00653400000000000	1.64956100000000	0.150362000000000	0.537072000000000	0.0319140000000000	8.00995800000000	201.800293000000	24.0394820000000
    'Ti3+'	22	36.5879330000000	7.23025500000000	-9.08607700000000	2.08459400000000	17.2940080000000	-35.1112820000000	0.000681000000000000	0.522262000000000	5.26231700000000	15.8817160000000	6.14980500000000
    'Ti4+'	22	45.3555370000000	7.09290000000000	7.48385800000000	-43.4988170000000	1.67891500000000	-0.110628000000000	9.25218600000000	0.523046000000000	13.0828520000000	10.1938760000000	0.0230640000000000
    'V'	    23	10.4735750000000	1.54788100000000	1.98638100000000	1.86561600000000	7.05625000000000	0.0677440000000000	7.08194000000000	0.0260400000000000	31.9096720000000	108.022842000000	0.474882000000000
    'V2+'	23	7.75435600000000	2.06410000000000	2.57699800000000	2.01140400000000	7.12617700000000	-0.533379000000000	7.06631500000000	0.0149930000000000	7.06630800000000	22.0557860000000	0.467568000000000
    'V3+'	23	9.95848000000000	1.59635000000000	1.48344200000000	-10.8460440000000	17.3328670000000	0.474921000000000	6.76304100000000	0.0568950000000000	17.7500290000000	0.328826000000000	0.388013000000000
    'V5+'	23	15.5750180000000	8.44809500000000	1.61204000000000	-9.72185500000000	1.53402900000000	0.552676000000000	0.682708000000000	5.56664000000000	10.5270770000000	0.907961000000000	0.0666670000000000
    'Cr'	24	11.0070690000000	1.55547700000000	2.98529300000000	1.34785500000000	7.03477900000000	0.0655100000000000	6.36628100000000	0.0239870000000000	23.2448390000000	105.774498000000	0.429369000000000
    'Cr2+'	24	10.5988770000000	1.56585800000000	2.72828000000000	0.0980640000000000	6.95932100000000	0.0498700000000000	6.15184600000000	0.0235190000000000	17.4328160000000	54.0023880000000	0.426301000000000
    'Cr3+'	24	7.98931000000000	1.76507900000000	2.62712500000000	1.82938000000000	6.98090800000000	-0.192123000000000	6.06886700000000	0.0183420000000000	6.06888700000000	16.3092840000000	0.420864000000000
    'Mn'	25	11.7095420000000	1.73341400000000	2.67314100000000	2.02336800000000	7.00318000000000	-0.147293000000000	5.59712000000000	0.0178000000000000	21.7884200000000	89.5179140000000	0.383054000000000
    'Mn2+'	25	11.2877120000000	26.0424140000000	3.05809600000000	0.0902580000000000	7.08830600000000	-24.5661320000000	5.50622500000000	0.000774000000000000	16.1585750000000	54.7663540000000	0.375580000000000
    'Mn3+'	25	6.92697200000000	2.08134200000000	11.1283790000000	2.37510700000000	-0.419287000000000	-0.0937130000000000	0.378315000000000	0.0150540000000000	5.37995700000000	14.4295860000000	0.00493900000000000
    'Mn4+'	25	12.4091310000000	7.46699300000000	1.80994700000000	-12.1384770000000	10.7802480000000	0.672146000000000	0.300400000000000	0.112814000000000	12.5207560000000	0.168653000000000	5.17323700000000
    'Fe'	26	12.3110980000000	1.87662300000000	3.06617700000000	2.07045100000000	6.97518500000000	-0.304931000000000	5.00941500000000	0.0144610000000000	18.7430400000000	82.7678760000000	0.346506000000000
    'Fe2+'	26	11.7767650000000	11.1650970000000	3.53349500000000	0.165345000000000	7.03693200000000	-9.67691900000000	4.91223200000000	0.00174800000000000	14.1665560000000	42.3819580000000	0.341324000000000
    'Fe3+'	26	9.72163800000000	63.4038470000000	2.14134700000000	2.62927400000000	7.03384600000000	-61.9307250000000	4.86929700000000	0.000293000000000000	4.86760200000000	13.5390760000000	0.338520000000000
    'Co'	27	12.9145100000000	2.48190800000000	3.46689400000000	2.10635100000000	6.96089200000000	-0.936572000000000	4.50713800000000	0.00912600000000000	16.4381290000000	76.9873200000000	0.314418000000000
    'Co2+'	27	6.99384000000000	26.2858120000000	12.2542890000000	0.246114000000000	4.01740700000000	-24.7968520000000	0.310779000000000	0.000684000000000000	4.40052800000000	35.7414470000000	12.5363930000000
    'Co3+'	27	6.86173900000000	2.67857000000000	12.2818890000000	3.50174100000000	-0.179384000000000	-1.14734500000000	0.309794000000000	0.00814200000000000	4.33170300000000	11.9141670000000	11.9141670000000
    'Ni'	28	13.5218650000000	6.94728500000000	3.86602800000000	2.13590000000000	4.28473100000000	-2.76269700000000	4.07727700000000	0.286763000000000	14.6226340000000	71.9660800000000	0.00443700000000000
    'Ni2+'	28	12.5190170000000	37.8320580000000	4.38725700000000	0.661552000000000	6.94907200000000	-36.3444710000000	3.93305300000000	0.000442000000000000	10.4491840000000	23.8609980000000	0.283723000000000
    'Ni3+'	28	13.5793660000000	1.90284400000000	12.8592680000000	3.81100500000000	-6.83859500000000	-0.317618000000000	0.313140000000000	0.0126210000000000	3.90640700000000	10.8943110000000	0.344379000000000
    'Cu'	29	14.0141920000000	4.78457700000000	5.05680600000000	1.45797100000000	6.93299600000000	-3.25447700000000	3.73828000000000	0.00374400000000000	13.0349820000000	72.5547940000000	0.265666000000000
    'Cu1+'	29	12.9607630000000	16.3421500000000	1.11010200000000	5.52068200000000	6.91545200000000	-14.8493200000000	3.57601000000000	0.000975000000000000	29.5232180000000	10.1142830000000	0.261326000000000
    'Cu2+'	29	11.8955690000000	16.3449780000000	5.79981700000000	1.04880400000000	6.78908800000000	-14.8783830000000	3.37851900000000	0.000924000000000000	8.13365300000000	20.5265240000000	0.254741000000000
    'Zn'	30	14.7410020000000	6.90774800000000	4.64233700000000	2.19176600000000	38.4240420000000	-36.9158290000000	3.38823200000000	0.243315000000000	11.9036890000000	63.3121300000000	0.000397000000000000
    'Zn2+'	30	13.3407720000000	10.4288570000000	5.54448900000000	0.762295000000000	6.86917200000000	-8.94524800000000	3.21591300000000	0.00141300000000000	8.54268000000000	21.8917560000000	0.239215000000000
    'Ga'	31	15.7589460000000	6.84112300000000	4.12101600000000	2.71468100000000	2.39524600000000	-0.847395000000000	3.12175400000000	0.226057000000000	12.4821960000000	66.2036210000000	0.00723800000000000
    'Ga3+'	31	13.1238750000000	35.2881890000000	6.12697900000000	0.611551000000000	6.72480700000000	-33.8751220000000	2.80996000000000	0.000323000000000000	6.83153400000000	16.7843110000000	0.212002000000000
    'Ge'	32	16.5406130000000	1.56790000000000	3.72782900000000	3.34509800000000	6.78507900000000	0.0187260000000000	2.86661800000000	0.0121980000000000	13.4321630000000	58.8660470000000	0.210974000000000
    'Ge4+'	32	6.87663600000000	6.77909100000000	9.96959100000000	3.13585700000000	0.152389000000000	1.08654200000000	2.02517400000000	0.176650000000000	3.57382200000000	7.68584800000000	16.6775740000000
    'As'	33	17.0256420000000	4.50344100000000	3.71590400000000	3.93720000000000	6.79017500000000	-2.98411700000000	2.59773900000000	0.00301200000000000	14.2721190000000	50.4379960000000	0.193015000000000
    'Se'	34	17.3540710000000	4.65324800000000	4.25948900000000	4.13645500000000	6.74916300000000	-3.16098200000000	2.34978700000000	0.00255000000000000	15.5794600000000	45.1812020000000	0.177432000000000
    'Br'	35	17.5505700000000	5.41188200000000	3.93718000000000	3.88064500000000	6.70779300000000	-2.49208800000000	2.11922600000000	16.5571840000000	0.00248100000000000	42.1640090000000	0.162121000000000
    'Br1-'	35	17.7143100000000	6.46692600000000	6.94738500000000	4.40267400000000	-0.697279000000000	1.15267400000000	2.12255400000000	19.0507680000000	0.152708000000000	58.6903610000000	58.6903720000000
    'Kr'	36	17.6552790000000	6.84810500000000	4.17100400000000	3.44676000000000	6.68520000000000	-2.81059200000000	1.90823100000000	16.6062360000000	0.00159800000000000	39.9174730000000	0.146896000000000
    'Rb'	37	8.12313400000000	2.13804200000000	6.76170200000000	1.15605100000000	17.6795460000000	1.13954800000000	15.1423850000000	33.5426670000000	0.129372000000000	224.132507000000	1.71336800000000
    'Rb1+'	37	17.6843200000000	7.76158800000000	6.68087400000000	2.66888300000000	0.0709740000000000	1.13326300000000	1.71020900000000	14.9198630000000	0.128542000000000	31.6544780000000	0.128543000000000
    'Sr'	38	17.7302190000000	9.79586700000000	6.09976300000000	2.62002500000000	0.600053000000000	1.14025100000000	1.56306000000000	14.3108680000000	0.120574000000000	135.771317000000	0.120574000000000
    'Sr2+'	38	17.6949730000000	1.27576200000000	6.15425200000000	9.23478600000000	0.515995000000000	1.12530900000000	1.55088800000000	30.1330410000000	0.118774000000000	13.8217990000000	0.118774000000000
    'Y'	    39	17.7920400000000	10.2532520000000	5.71494900000000	3.17051600000000	0.918251000000000	1.13178700000000	1.42969100000000	13.1328160000000	0.112173000000000	108.197029000000	0.112173000000000
    'Zr'	40	17.8597720000000	10.9110380000000	5.82111500000000	3.51251300000000	0.746965000000000	1.12485900000000	1.31069200000000	12.3192850000000	0.104353000000000	91.7775420000000	0.104353000000000
    'Zr4+'	40	6.80295600000000	17.6992530000000	10.6506470000000	-0.248108000000000	0.250338000000000	0.827902000000000	0.0962280000000000	1.29612700000000	11.2407150000000	-0.219259000000000	-0.219021000000000
    'Nb'	41	17.9583990000000	12.0630540000000	5.00701500000000	3.28766700000000	1.53101900000000	1.12345200000000	1.21159000000000	12.2466870000000	0.0986150000000000	75.0119480000000	0.0986150000000000
    'Nb3+'	41	17.7143230000000	1.67521300000000	7.48396300000000	8.32246400000000	11.1435730000000	-8.33957300000000	1.17241900000000	30.1027910000000	0.0802550000000000	-0.00298300000000000	10.4566870000000
    'Nb5+'	41	17.5802060000000	7.63327700000000	10.7934970000000	0.180884000000000	67.8379210000000	-68.0247800000000	1.16585200000000	0.0785580000000000	9.50765200000000	31.6216560000000	-0.000438000000000000
    'Mo'	42	6.23621800000000	17.9877110000000	12.9731270000000	3.45142600000000	0.210899000000000	1.10877000000000	0.0907800000000000	1.10831000000000	11.4687200000000	66.6841510000000	0.0907800000000000
    'Mo3+'	42	7.44705000000000	17.7781220000000	11.8860680000000	1.99790500000000	1.78962600000000	-1.89876400000000	0.0720000000000000	1.07314500000000	9.83472000000000	28.2217460000000	-0.0116740000000000
    'Mo5+'	42	7.92987900000000	17.6676690000000	11.5159870000000	0.500402000000000	77.4440840000000	-78.0565950000000	0.0688560000000000	1.06806400000000	9.04622900000000	26.5589450000000	-0.000473000000000000
    'Mo6+'	42	34.7576830000000	9.65303700000000	6.58476900000000	-18.6281150000000	2.49059400000000	1.14191600000000	1.30177000000000	7.12384300000000	0.0940970000000000	1.61744300000000	12.3354340000000
    'Tc'	43	17.8409630000000	3.42823600000000	1.37301200000000	12.9473640000000	6.33546900000000	1.07478400000000	1.00572900000000	41.9013820000000	119.320541000000	9.78154200000000	0.0833910000000000
    'Ru'	44	6.27162400000000	17.9067380000000	14.1232690000000	3.74600800000000	0.908235000000000	1.04399200000000	0.0770400000000000	0.928222000000000	9.55534500000000	35.8606800000000	123.552246000000
    'Ru3+'	44	17.8947580000000	13.5795290000000	10.7292510000000	2.47409500000000	48.2279970000000	-51.9052430000000	0.902827000000000	8.74057900000000	0.0451250000000000	24.7649540000000	-0.00169900000000000
    'Ru4+'	44	17.8457760000000	13.4550840000000	10.2290870000000	1.65352400000000	14.0597950000000	-17.2417620000000	0.901070000000000	8.48239200000000	0.0459720000000000	23.0152720000000	-0.00488900000000000
    'Rh'	45	6.21664800000000	17.9197390000000	3.85425200000000	0.840326000000000	15.1734980000000	0.995452000000000	0.0707890000000000	0.856121000000000	33.8894840000000	121.686691000000	9.02951700000000
    'Rh3+'	45	17.7586210000000	14.5698130000000	5.29832000000000	2.53357900000000	0.879753000000000	0.960843000000000	0.841779000000000	8.31953300000000	0.0690500000000000	23.7091310000000	0.0690500000000000
    'Rh4+'	45	17.7161880000000	14.4466540000000	5.18580100000000	1.70344800000000	0.989992000000000	0.959941000000000	0.840572000000000	8.10064700000000	0.0689950000000000	22.3573070000000	0.0689950000000000
    'Pd'	46	6.12151100000000	4.78406300000000	16.6316830000000	4.31825800000000	13.2467730000000	0.883099000000000	0.0625490000000000	0.784031000000000	8.75139100000000	34.4899830000000	0.784031000000000
    'Pd2+'	46	6.12228200000000	15.6510120000000	3.51350800000000	9.06079000000000	8.77119900000000	0.879336000000000	0.0624240000000000	8.01829600000000	24.7842750000000	0.776457000000000	0.776457000000000
    'Pd4+'	46	6.15242100000000	-96.0690230000000	31.6221410000000	81.5782550000000	17.8014030000000	0.915874000000000	0.0639510000000000	11.0903540000000	13.4661520000000	9.75830200000000	0.783014000000000
    'Ag'	47	6.07387400000000	17.1554370000000	4.17334400000000	0.852238000000000	17.9886860000000	0.756603000000000	0.0553330000000000	7.89651200000000	28.4437390000000	110.376106000000	0.716809000000000
    'Ag1+'	47	6.09119200000000	4.01952600000000	16.9481740000000	4.25863800000000	13.8894370000000	0.785127000000000	0.0563050000000000	0.719340000000000	7.75893800000000	27.3683490000000	0.719340000000000
    'Ag2+'	47	6.40180800000000	48.6998020000000	4.79985900000000	-32.3325230000000	16.3567100000000	1.06824700000000	0.0681670000000000	0.942270000000000	20.6394960000000	1.10036500000000	6.88313100000000
    'Cd'	48	6.08098600000000	18.0194680000000	4.01819700000000	1.30351000000000	17.9746690000000	0.603504000000000	0.0489900000000000	7.27364600000000	29.1192840000000	95.8312070000000	0.661231000000000
    'Cd2+'	48	6.09371100000000	43.9096910000000	17.0413060000000	-39.6751170000000	17.9589180000000	0.664795000000000	0.0506240000000000	8.65414300000000	15.6213960000000	11.0820670000000	0.667591000000000
    'In'	49	6.19647700000000	18.8161830000000	4.05047900000000	1.63892900000000	17.9629120000000	0.333097000000000	0.0420720000000000	6.69566500000000	31.0097900000000	103.284348000000	0.610714000000000
    'In3+'	49	6.20627700000000	18.4977460000000	3.07813100000000	10.5246130000000	7.40123400000000	0.293677000000000	0.0413570000000000	6.60556300000000	18.7922500000000	0.608082000000000	0.608082000000000
    'Sn'	50	19.3251710000000	6.28157100000000	4.49886600000000	1.85693400000000	17.9173180000000	0.119024000000000	6.11810400000000	0.0369150000000000	32.5290450000000	95.0371860000000	0.565651000000000
    'Sn2+'	50	6.35367200000000	4.77037700000000	14.6720250000000	4.23595900000000	18.0021310000000	-0.0425190000000000	0.0347200000000000	6.16789100000000	6.16787900000000	29.0064560000000	0.561774000000000
    'Sn4+'	50	15.4457320000000	6.42089200000000	4.56298000000000	1.71338500000000	18.0335370000000	-0.172219000000000	6.28089800000000	0.0331440000000000	6.28089900000000	17.9836010000000	0.557980000000000
    'Sb'	51	5.39495600000000	6.54957000000000	19.6506810000000	1.82782000000000	17.8678320000000	-0.290506000000000	33.3265230000000	0.0309740000000000	5.56492900000000	87.1309660000000	0.523992000000000
    'Sb3+'	51	10.1891710000000	57.4619180000000	19.3565730000000	4.86220600000000	-45.3940960000000	1.51610800000000	0.0894850000000000	0.375256000000000	5.35798700000000	22.1537360000000	0.297768000000000
    'Sb5+'	51	17.9206220000000	6.64793200000000	12.7240750000000	1.55554500000000	7.60059100000000	-0.445371000000000	0.522315000000000	0.0294870000000000	5.71821000000000	16.4337750000000	5.71820400000000
    'Te'	52	6.66030200000000	6.94075600000000	19.8470150000000	1.55717500000000	17.8024270000000	-0.806668000000000	33.0316540000000	0.0257500000000000	5.06554700000000	84.1016160000000	0.487660000000000
    'I' 	53	19.8845020000000	6.73659300000000	8.11051600000000	1.17095300000000	17.5487160000000	-0.448811000000000	4.62859100000000	0.0277540000000000	31.8490960000000	84.4063870000000	0.463550000000000
    'I1-'	53	20.0103300000000	17.8355240000000	8.10413000000000	2.23111800000000	9.15854800000000	-3.34100400000000	4.56593100000000	0.444266000000000	32.4306720000000	95.1490400000000	0.0149060000000000
    'Xe'	54	19.9789200000000	11.7749450000000	9.33218200000000	1.24474900000000	17.7375010000000	-6.06590200000000	4.14335600000000	0.0101420000000000	28.7962000000000	75.2806850000000	0.413616000000000
    'Cs'	55	17.4186740000000	8.31444400000000	10.3231930000000	1.38383400000000	19.8762510000000	-2.32280200000000	0.399828000000000	0.0168720000000000	25.6058270000000	233.339676000000	3.82691500000000
    'Cs1+'	55	19.9390560000000	24.9676210000000	10.3758840000000	0.454243000000000	17.6602480000000	-19.3943060000000	3.77051100000000	0.00404000000000000	25.3112750000000	76.5377660000000	0.384730000000000
    'Ba'	56	19.7473430000000	17.3684770000000	10.4657180000000	2.59260200000000	11.0036530000000	-5.18349700000000	3.48182300000000	0.371224000000000	21.2266410000000	173.834274000000	0.0107190000000000
    'Ba2+'	56	19.7502000000000	17.5136830000000	10.8848920000000	0.321585000000000	65.1498340000000	-59.6181720000000	3.43074800000000	0.361590000000000	21.3583070000000	70.3094020000000	0.00141800000000000
    'La'	57	19.9660190000000	27.3296550000000	11.0184250000000	3.08669600000000	17.3354550000000	-21.7454890000000	3.19740800000000	0.00344600000000000	19.9554920000000	141.381973000000	0.341817000000000
    'La3+'	57	19.6888870000000	17.3457030000000	11.3562960000000	0.0994180000000000	82.3581240000000	-76.8469090000000	3.14621100000000	0.339586000000000	18.7538320000000	90.3454590000000	0.00107200000000000
    'Ce'	58	17.3551220000000	43.9884990000000	20.5466500000000	3.13067000000000	11.3536650000000	-38.3860170000000	0.328369000000000	0.00204700000000000	3.08819600000000	134.907654000000	18.8329600000000
    'Ce3+'	58	26.5932310000000	85.8664320000000	-6.67769500000000	12.1118470000000	17.4019030000000	-80.3134230000000	3.28038100000000	0.00101200000000000	4.31357500000000	17.8685040000000	0.326962000000000
    'Ce4+'	58	17.4575330000000	25.6599410000000	11.6910370000000	19.6952510000000	-16.9947490000000	-3.51509600000000	0.311812000000000	-0.00379300000000000	16.5686870000000	2.88639500000000	-0.00893100000000000
    'Pr'	59	21.5513110000000	17.1617300000000	11.9038590000000	2.67910300000000	9.56419700000000	-3.87106800000000	2.99567500000000	0.312491000000000	17.7167050000000	152.192825000000	0.0104680000000000
    'Pr3+'	59	20.8798410000000	36.0357970000000	12.1353410000000	0.283103000000000	17.1678030000000	-30.5007840000000	2.87089700000000	0.00236400000000000	16.6152360000000	53.9093590000000	0.306993000000000
    'Pr4+'	59	17.4960820000000	21.5385090000000	20.4031140000000	12.0622110000000	-7.49204300000000	-9.01672200000000	0.294457000000000	-0.00274200000000000	2.77288600000000	15.8046130000000	-0.0135560000000000
    'Nd'	60	17.3312440000000	62.7839240000000	12.1600970000000	2.66348300000000	22.2399500000000	-57.1898420000000	0.300269000000000	0.00132000000000000	17.0260010000000	148.748993000000	2.91026800000000
    'Nd3+'	60	17.1200770000000	56.0381390000000	21.4683070000000	10.0006710000000	2.90586600000000	-50.5419920000000	0.291295000000000	0.00142100000000000	2.74368100000000	14.5813670000000	22.4850980000000
    'Pm'	61	17.2863880000000	51.5601620000000	12.4785570000000	2.67551500000000	22.9609470000000	-45.9736820000000	0.286620000000000	0.00155000000000000	16.2237550000000	143.984512000000	2.79648000000000
    'Pm3+'	61	22.2210660000000	17.0681420000000	12.8054230000000	0.435687000000000	52.2387700000000	-46.7671810000000	2.63576700000000	0.277039000000000	14.9273150000000	45.7680170000000	0.00145500000000000
    'Sm'	62	23.7003630000000	23.0722140000000	12.7777820000000	2.68421700000000	17.2043670000000	-17.4521660000000	2.68953900000000	0.00349100000000000	15.4954370000000	139.862473000000	0.274536000000000
    'Sm3+'	62	15.6185650000000	19.5380920000000	13.3989460000000	-4.35881100000000	24.4904610000000	-9.71485400000000	0.00600100000000000	0.306379000000000	14.9795940000000	0.748825000000000	2.45449200000000
    'Eu'	63	17.1861950000000	37.1568370000000	13.1033870000000	2.70724600000000	24.4192710000000	-31.5866870000000	0.261678000000000	0.00199500000000000	14.7873600000000	134.816299000000	2.58188300000000
    'Eu2+'	63	23.8990350000000	31.6574970000000	12.9557520000000	1.70057600000000	16.9921990000000	-26.2043150000000	2.46733200000000	0.00223000000000000	13.6250020000000	35.0894810000000	0.253136000000000
    'Eu3+'	63	17.7583270000000	33.4986650000000	24.0671880000000	13.4368830000000	-9.01913400000000	-19.7680260000000	0.244474000000000	-0.00390100000000000	2.48752600000000	14.5680110000000	-0.0156280000000000
    'Gd'	64	24.8981170000000	17.1049520000000	13.2225810000000	3.26615200000000	48.9952130000000	-43.5056840000000	2.43502800000000	0.246961000000000	13.9963250000000	110.863091000000	0.00138300000000000
    'Gd3+'	64	24.3449990000000	16.9453110000000	13.8669310000000	0.481674000000000	93.5063780000000	-88.1471790000000	2.33397100000000	0.239215000000000	12.9829950000000	43.8763470000000	0.000673000000000000
    'Tb'	65	25.9100130000000	32.3441390000000	13.7651170000000	2.75140400000000	17.0644050000000	-26.8519710000000	2.37391200000000	0.00203400000000000	13.4819690000000	125.836510000000	0.236916000000000
    'Tb3+'	65	24.8782520000000	16.8560160000000	13.6639370000000	1.27967100000000	39.2712940000000	-33.9503170000000	2.22330100000000	0.227290000000000	11.8125280000000	29.9100650000000	0.00152700000000000
    'Dy'	66	26.6717850000000	88.6875760000000	14.0654450000000	2.76849700000000	17.0677810000000	-83.2798310000000	2.28259300000000	0.000665000000000000	12.9202300000000	121.937187000000	0.225531000000000
    'Dy3+'	66	16.8643440000000	90.3834610000000	13.6754730000000	1.68707800000000	25.5406510000000	-85.1506500000000	0.216275000000000	0.000593000000000000	11.1212070000000	26.2509750000000	2.13593000000000
    'Ho'	67	27.1501900000000	16.9998190000000	14.0593340000000	3.38697900000000	46.5464710000000	-41.1652530000000	2.16966000000000	0.215414000000000	12.2131480000000	100.506783000000	0.00121100000000000
    'Ho3+'	67	16.8375240000000	63.2213360000000	13.7037660000000	2.06160200000000	26.2026210000000	-58.0265050000000	0.206873000000000	0.000796000000000000	10.5002830000000	24.0318830000000	2.05506000000000
    'Er'	68	28.1748870000000	82.4932710000000	14.6240020000000	2.80275600000000	17.0185150000000	-77.1352230000000	2.12099500000000	0.000640000000000000	11.9152560000000	114.529938000000	0.207519000000000
    'Er3+'	68	16.8101270000000	22.6810610000000	13.8641140000000	2.29450600000000	26.8644770000000	-17.5134600000000	0.198293000000000	0.00212600000000000	9.97334100000000	22.8363880000000	1.97944200000000
    'Tm'	69	28.9258940000000	76.1737980000000	14.9047040000000	2.81481200000000	16.9981170000000	-70.8398130000000	2.04620300000000	0.000656000000000000	11.4653750000000	111.411980000000	0.199376000000000
    'Tm3+'	69	16.7875000000000	15.3509050000000	14.1823570000000	2.29911100000000	27.5737710000000	-10.1920870000000	0.190852000000000	0.00303600000000000	9.60293400000000	22.5268800000000	1.91286200000000
    'Yb'	70	29.6767600000000	65.6240690000000	15.1608540000000	2.83028800000000	16.9978500000000	-60.3138120000000	1.97763000000000	0.000720000000000000	11.0446220000000	108.139153000000	0.192110000000000
    'Yb2+'	70	28.4437940000000	16.8495270000000	14.1650810000000	3.44531100000000	28.3088530000000	-23.2149350000000	1.86389600000000	0.183811000000000	9.22546900000000	23.6913550000000	0.00146300000000000
    'Yb3+'	70	28.1916290000000	16.8280870000000	14.1678480000000	2.74496200000000	23.1717740000000	-18.1036760000000	1.84288900000000	0.182788000000000	9.04595700000000	20.7998470000000	0.00175900000000000
    'Lu'	71	30.1228660000000	15.0993460000000	56.3148990000000	3.54098000000000	16.9437290000000	-51.0494160000000	1.88309000000000	10.3427640000000	0.000780000000000000	89.5592500000000	0.183849000000000
    'Lu3+'	71	28.8286930000000	16.8232270000000	14.2476170000000	3.07955900000000	25.6476670000000	-20.6265280000000	1.77664100000000	0.175560000000000	8.57553100000000	19.6937010000000	0.00145300000000000
    'Hf'	72	30.6170330000000	15.1453510000000	54.9335480000000	4.09625300000000	16.8961560000000	-49.7198370000000	1.79561300000000	9.93446900000000	0.000739000000000000	76.1897050000000	0.175914000000000
    'Hf4+'	72	29.2673780000000	16.7925430000000	14.7853100000000	2.18412800000000	23.7919960000000	-18.8203830000000	1.69791100000000	0.168313000000000	8.19002500000000	18.2775780000000	0.00143100000000000
    'Ta'	73	31.0663590000000	15.3418230000000	49.2782970000000	4.57766500000000	16.8283210000000	-44.1190260000000	1.70873200000000	9.61845500000000	0.000760000000000000	66.3461990000000	0.168002000000000
    'Ta5+'	73	29.5394690000000	16.7418540000000	15.1820700000000	1.64291600000000	16.4374470000000	-11.5424590000000	1.61293400000000	0.160460000000000	7.65440800000000	17.0707320000000	0.00185800000000000
    'W' 	74	31.5079000000000	15.6824980000000	37.9601290000000	4.88550900000000	16.7921120000000	-32.8645740000000	1.62948500000000	9.44644800000000	0.000898000000000000	59.9806750000000	0.160798000000000
    'W6+'	74	29.7293570000000	17.2478080000000	15.1844880000000	1.15465200000000	0.739335000000000	3.94515700000000	1.50164800000000	0.140803000000000	6.88057300000000	14.2996010000000	14.2996180000000
    'Re'	75	31.8884560000000	16.1171040000000	42.3902970000000	5.21166900000000	16.7675910000000	-37.4126820000000	1.54923800000000	9.23347400000000	0.000689000000000000	54.5163730000000	0.152815000000000
    'Os'	76	32.2102970000000	16.6784400000000	48.5599060000000	5.45583900000000	16.7355330000000	-43.6779560000000	1.47353100000000	9.04969500000000	0.000519000000000000	50.2102010000000	0.145771000000000
    'Os4+'	76	17.1134850000000	15.7923700000000	23.3423920000000	4.09027100000000	7.67129200000000	3.98839000000000	0.131850000000000	7.28854200000000	1.38930700000000	19.6294250000000	1.38930700000000
    'Ir'	77	32.0044360000000	1.97545400000000	17.0701050000000	15.9394540000000	5.99000300000000	4.01889300000000	1.35376700000000	81.0141750000000	0.128093000000000	7.66119600000000	26.6594030000000
    'Ir3+'	77	31.5375750000000	16.3633380000000	15.5971410000000	5.05140400000000	1.43693500000000	4.00945900000000	1.33414400000000	7.45191800000000	0.127514000000000	21.7056480000000	0.127515000000000
    'Ir4+'	77	30.3912490000000	16.1469960000000	17.0190680000000	4.45890400000000	0.975372000000000	4.00686500000000	1.32851900000000	7.18176600000000	0.127337000000000	19.0601460000000	1.32851900000000
    'Pt'	78	31.2738910000000	18.4454400000000	17.0637450000000	5.55593300000000	1.57527000000000	4.05039400000000	1.31699200000000	8.79715400000000	0.124741000000000	40.1779940000000	1.31699700000000
    'Pt2+'	78	31.9868490000000	17.2490480000000	15.2693740000000	5.76023400000000	1.69407900000000	4.03251200000000	1.28114300000000	7.62551200000000	0.123571000000000	24.1908260000000	0.123571000000000
    'Pt4+'	78	41.9327130000000	16.3392240000000	17.6538940000000	6.01242000000000	-12.0368770000000	4.09455100000000	1.11140900000000	6.46608600000000	0.128917000000000	16.9541550000000	0.778721000000000
    'Au'	79	16.7773900000000	19.3171560000000	32.9796830000000	5.59545300000000	10.5768540000000	-6.27907800000000	0.122737000000000	8.62157000000000	1.25690200000000	38.0088200000000	0.000601000000000000
    'Au1+'	79	32.1243060000000	16.7164760000000	16.8141000000000	7.31156500000000	0.993064000000000	4.04079200000000	1.21607300000000	7.16537800000000	0.118715000000000	20.4424860000000	53.0959850000000
    'Au3+'	79	31.7042710000000	17.5457670000000	16.8195510000000	5.52264000000000	0.361725000000000	4.04267900000000	1.21556100000000	7.22050600000000	0.118812000000000	20.0509700000000	1.21556200000000
    'Hg'	80	16.8398900000000	20.0238230000000	28.4285640000000	5.88156400000000	4.71470600000000	4.07647800000000	0.115905000000000	8.25692700000000	1.19525000000000	39.2472270000000	1.19525000000000
    'Hg1+'	80	28.8668370000000	19.2775400000000	16.7760510000000	6.28145900000000	3.71028900000000	4.06843000000000	1.17396700000000	7.58384200000000	0.115351000000000	29.0559940000000	1.17396800000000
    'Hg2+'	80	32.4110790000000	18.6903710000000	16.7117730000000	9.97483500000000	-3.84761100000000	4.05286900000000	1.16298000000000	7.32980600000000	0.114518000000000	22.0094890000000	22.0094930000000
    'Tl'	81	16.6307950000000	19.3866160000000	32.8085710000000	1.74719100000000	6.35686200000000	4.06693900000000	0.110704000000000	7.18140100000000	1.11973000000000	90.6602630000000	26.0149780000000
    'Tl1+'	81	32.2950440000000	16.5700490000000	17.9910130000000	1.53535500000000	7.55459100000000	4.05403000000000	1.10154400000000	0.110020000000000	6.52855900000000	52.4950680000000	20.3386340000000
    'Tl3+'	81	32.5256390000000	19.1391850000000	17.1003210000000	5.89111500000000	12.5994630000000	-9.25607500000000	1.09496600000000	6.90099200000000	0.103667000000000	18.4896140000000	-0.00140100000000000
    'Pb'	82	16.4195670000000	32.7385900000000	6.53024700000000	2.34274200000000	19.9164750000000	4.04982400000000	0.105499000000000	1.05504900000000	25.0258900000000	80.9065930000000	6.66444900000000
    'Pb2+'	82	27.3926470000000	16.4968220000000	19.9845010000000	6.81392300000000	5.23391000000000	4.06562300000000	1.05887400000000	0.106305000000000	6.70812300000000	24.3955540000000	1.05887400000000
    'Pb4+'	82	32.5056570000000	20.0142400000000	14.6456610000000	5.02949900000000	1.76013800000000	4.04467800000000	1.04703500000000	6.67032100000000	0.105279000000000	16.5250400000000	0.105279000000000
    'Bi'	83	16.2822740000000	32.7251360000000	6.67830200000000	2.69475000000000	20.5765590000000	4.04091400000000	0.101180000000000	1.00228700000000	25.7141460000000	77.0575490000000	6.29188200000000
    'Bi3+'	83	32.4614370000000	19.4386830000000	16.3024860000000	7.32266200000000	0.431704000000000	4.04370300000000	0.997930000000000	6.03886700000000	0.101338000000000	18.3715860000000	46.3610460000000
    'Bi5+'	83	16.7340280000000	20.5804940000000	9.45262300000000	61.1558340000000	-34.0410230000000	4.11366300000000	0.105076000000000	4.77328200000000	11.7621620000000	1.21177500000000	1.61940800000000
    'Po'	84	16.2891640000000	32.8071710000000	21.0951630000000	2.50590100000000	7.25458900000000	4.04655600000000	0.0981210000000000	0.966265000000000	6.04662200000000	76.5980680000000	28.0961280000000
    'At'	85	16.0114610000000	32.6155470000000	8.11389900000000	2.88408200000000	21.3778670000000	3.99568400000000	0.0926390000000000	0.904416000000000	26.5432570000000	68.3729630000000	5.49951200000000
    'Rn'	86	16.0702290000000	32.6411060000000	21.4896580000000	2.29921800000000	9.48018400000000	4.02097700000000	0.0904370000000000	0.876409000000000	5.23968700000000	69.1884770000000	27.6326410000000
    'Fr'	87	16.0073850000000	32.6638300000000	21.5943510000000	1.59849700000000	11.1211920000000	4.00347200000000	0.0870310000000000	0.840187000000000	4.95446700000000	199.805801000000	26.9051060000000
    'Ra'	88	32.5636900000000	21.3966710000000	11.2980930000000	2.83468800000000	15.9149650000000	3.98177300000000	0.801980000000000	4.59066600000000	22.7589720000000	160.404388000000	0.0835440000000000
    'Ra2+'	88	4.98622800000000	32.4749450000000	21.9474430000000	11.8000130000000	10.8072920000000	3.95657200000000	0.0825970000000000	0.791468000000000	4.60803400000000	24.7924310000000	0.0825970000000000
    'Ac'	89	15.9140530000000	32.5350420000000	21.5539760000000	11.4333940000000	3.61240900000000	3.93921200000000	0.0805110000000000	0.770669000000000	4.35220600000000	21.3816220000000	130.500748000000
    'Ac3+'	89	15.5849830000000	32.0221250000000	21.4563270000000	0.757593000000000	12.3412520000000	3.83898400000000	0.0774380000000000	0.739963000000000	4.04073500000000	47.5250020000000	19.4068450000000
    'Th'	90	15.7840240000000	32.4548990000000	21.8492220000000	4.23907700000000	11.7361910000000	3.92253300000000	0.0770670000000000	0.735137000000000	4.09797600000000	109.464111000000	20.5121380000000
    'Th4+'	90	15.5154450000000	32.0906910000000	13.9963990000000	12.9181570000000	7.63551400000000	3.83112200000000	0.0744990000000000	0.711663000000000	3.87104400000000	18.5968910000000	3.87104400000000
    'Pa'	91	32.7402080000000	21.9736750000000	12.9573980000000	3.68383200000000	15.7440580000000	3.88606600000000	0.709545000000000	4.05088100000000	19.2315430000000	117.255005000000	0.0740400000000000
    'U' 	92	15.6792750000000	32.8243060000000	13.6604590000000	3.68726100000000	22.2794340000000	3.85444400000000	0.0712060000000000	0.681177000000000	18.2361560000000	112.500038000000	3.93032500000000
    'U3+'	92	15.3603090000000	32.3956570000000	21.9612900000000	1.32589400000000	14.2514530000000	3.70662200000000	0.0678150000000000	0.654643000000000	3.64340900000000	39.6049650000000	16.3305700000000
    'U4+'	92	15.3550910000000	32.2353060000000	0.557745000000000	14.3963670000000	21.7511730000000	3.70586300000000	0.0677890000000000	0.652613000000000	42.3542370000000	15.9082390000000	3.55323100000000
    'U6+'	92	15.3338440000000	31.7708490000000	21.2744140000000	13.8726360000000	0.0485190000000000	3.70059100000000	0.0676440000000000	0.646384000000000	3.31789400000000	14.6502500000000	75.3396990000000
    'Np'	93	32.9999010000000	22.6380770000000	14.2199730000000	3.67295000000000	15.6832450000000	3.76939100000000	0.657086000000000	3.85491800000000	17.4354740000000	109.464485000000	0.0680330000000000
    'Np3+'	93	15.3781520000000	32.5721320000000	22.2061250000000	1.41329500000000	14.8283810000000	3.60337000000000	0.0646130000000000	0.631420000000000	3.56193600000000	37.8755110000000	15.5461290000000
    'Np4+'	93	15.3739260000000	32.4230190000000	21.9699940000000	0.662078000000000	14.9693500000000	3.60303900000000	0.0645970000000000	0.629658000000000	3.47638900000000	39.4389420000000	15.1357640000000
    'Np6+'	93	15.3599860000000	31.9928250000000	21.4124580000000	0.0665740000000000	14.5681740000000	3.60094200000000	0.0645280000000000	0.624505000000000	3.25344100000000	67.6583180000000	13.9808320000000
    'Pu'	94	33.2811780000000	23.1485440000000	15.1537550000000	3.03149200000000	15.7042150000000	3.66420000000000	0.634999000000000	3.85616800000000	16.8497350000000	121.292038000000	0.0648570000000000
    'Pu3+'	94	15.3560040000000	32.7691270000000	22.6802100000000	1.35105500000000	15.4162320000000	3.42889500000000	0.0605900000000000	0.604663000000000	3.49150900000000	37.2606350000000	14.9819210000000
    'Pu4+'	94	15.4162190000000	32.6105690000000	22.2566620000000	0.719495000000000	15.5181520000000	3.48040800000000	0.0614560000000000	0.607938000000000	3.41184800000000	37.6287920000000	14.4643600000000
    'Pu6+'	94	15.4365060000000	32.2897190000000	14.7267370000000	15.0123910000000	7.02467700000000	3.50232500000000	0.0618150000000000	0.606541000000000	3.24536300000000	13.6164380000000	3.24536400000000
    'Am'	95	33.4351620000000	23.6572590000000	15.5763390000000	3.02702300000000	15.7461000000000	3.54116000000000	0.612785000000000	3.79294200000000	16.1957780000000	117.757004000000	0.0617550000000000
    'Cm'	96	15.8048370000000	33.4808010000000	24.1501980000000	3.65556300000000	15.4998660000000	3.39084000000000	0.0586190000000000	0.590160000000000	3.67472000000000	100.736191000000	15.4082960000000
    'Bk'	97	15.8890720000000	33.6252860000000	24.7103810000000	3.70713900000000	15.8392680000000	3.21316900000000	0.0555030000000000	0.569571000000000	3.61547200000000	97.6947860000000	14.7543030000000
    'Cf'	98	33.7940750000000	25.4676930000000	16.0484870000000	3.65752500000000	16.0089820000000	3.00532600000000	0.550447000000000	3.58197300000000	14.3573880000000	96.0649720000000	0.0524500000000000 };

if isstring(name)
    name = char(name);
end
name = strtrim(name);

idx = find(strcmp(T(:,1), name), 1, 'first');

if isempty(idx)
    error('CMcoef:UnknownSpecies', 'CMcoef could not find coefficients for "%s".', name);
end

f = cell2mat(T(idx, 3:end)).';

if nargout > 1
    varargout{1} = T{idx, 2};
end
end



