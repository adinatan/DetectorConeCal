% RUN_LAB6AVG_TZ400_PILATUS3_2M
% LaB6 pose calibration for lab6avg_tz400.mat on a PILATUS3 X 2M.
%
% Uses:
%   fit_lab6_pilatus3_2m_oldlut_robust.m
%   cif2powder_CM.m
%   1000057.cif   (LaB6)
%
% Fitted parameter convention:
%   p = [beamX_px beamY_px zBeam_mm thetaX_deg thetaY_deg]
%
% Analysis of lab6avg_tz400.mat:
%   image size       = 1679 x 1475 [rows x cols]
%   beam intercept   ~ (49, 1605) px
%   distance prior   ~ 250-260 mm
%   detector tilt    ~ 0 deg
%
% IMPORTANT:
%   EkeV must be the actual X-ray energy for this image. 59.0 keV is used
%   here because that is the working value from the current calibration.
%   Replace it with the exact beamline value if different.

clear; clc; close all;
tTotal = tic;

%% ========================================================================
%  DATA AND BEAM ENERGY
% =========================================================================
dataFile = 'lab6avg_tz400.mat';
EkeV = 59.04;                        % <-- replace by exact value if different

%% ========================================================================
%  CONFIGURATION
% =========================================================================
cfg = struct();

% -------------------------------------------------------------------------
% PILATUS3 X 2M detector
% -------------------------------------------------------------------------
cfg.detector.name = 'PILATUS3_X_2M';
cfg.detector.pixel_um = 172;
cfg.detector.expectedSize_rc = [1679 1475];
cfg.detector.strictSize = true;

% -------------------------------------------------------------------------
% Beam
% -------------------------------------------------------------------------
cfg.beam.EkeV = EkeV;
cfg.beam.polMode = 'horizontal';

% lab6_avg is already 1679 rows x 1475 columns.
cfg.image.orientation = 'as_loaded';

% -------------------------------------------------------------------------
% Pose search bounds inferred from this image
%
% MATLAB image coordinates are 1-based:
%   x = column
%   y = row
%
% The q=0 intercept is near the lower-left corner at approximately
% (x,y) = (49,1605).
% -------------------------------------------------------------------------
cfg.bounds.beamX_px = [25 75];
cfg.bounds.beamY_px = [1575 1635];
cfg.bounds.zBeam_mm = [240 1080];
cfg.bounds.thetaX_deg = [-2 2];
cfg.bounds.thetaY_deg = [-2 2];

% Near-normal detector prior. This is only a soft guide; the bounds above
% remain the hard limits for the initial search.
cfg.prior.mainTiltAxis = 'y';
cfg.prior.tiltDeg = 0;
cfg.prior.tiltSign = +1;
cfg.prior.constrainBounds = false;
cfg.prior.tiltSigma_deg = 1.0;
cfg.prior.tiltPenaltyWeight = 0.002;
cfg.prior.penaltyType = 'soft';
cfg.prior.maxPenalty = 0.02;

% -------------------------------------------------------------------------
% LaB6 reference
%
% The old PILATUS100K code stopped near Q=9.5 A^-1. This PILATUS3 2M image
% spans much farther in Q, so explicitly extend the reference to 30 A^-1.
% -------------------------------------------------------------------------
cfg.ref.cifFile = '1000057.cif';
cfg.ref.altCifFile = '1000057(1).cif';
cfg.ref.qGrid_Ainv = 0.3:0.005:30.0;
cfg.ref.minRelPeakI = 0.004;
cfg.ref.profileFwhmQ_Ainv = 0.035;
cfg.ref.useCache = true;
cfg.ref.cacheFile = sprintf('lab6_reference_%.4gkeV_Q30_pilatus3_2m_cache.mat',EkeV);

% -------------------------------------------------------------------------
% Image preprocessing / ring-pixel selection
% -------------------------------------------------------------------------
cfg.prep.capHighPercentile = 99.98;
cfg.prep.backgroundWindowPx = 35;
cfg.prep.thresholdSigma = 2.5;
cfg.prep.thresholdPercentile = 98.5;
cfg.prep.minObjectPixels = 2;
cfg.prep.useClean = true;
cfg.prep.maxScorePixels = 6000;
cfg.prep.weightOffset = 0.25;
cfg.prep.computeLegacyAudit = false;

% Allow gentle threshold relaxation if the initial mask is too sparse.
cfg.prep.autoRelaxThreshold = true;
cfg.prep.minScorePixels = 1100;
cfg.prep.relaxSigmaFloor = 2.2;
cfg.prep.relaxPercentileFloor = 98.2;
cfg.prep.relaxSigmaStep = 0.15;
cfg.prep.relaxPercentileStep = 0.20;

% -------------------------------------------------------------------------
% CPU / GPU
% -------------------------------------------------------------------------
cfg.compute.device = 'auto';             % 'auto', 'cpu', or 'gpu'
cfg.compute.autoBenchmark = true;
cfg.compute.autoBenchmarkCandidates = 2048;
cfg.compute.autoGpuMinSpeedup = 1.05;
cfg.compute.cpuBlockCandidates = 512;
cfg.compute.gpuBlockCandidates = 4096;
cfg.compute.gpuPrecision = 'single';      % LUT only; final fit is double
cfg.compute.cpuParallel = 'auto';
cfg.compute.hitLookupDQ_Ainv = 1e-5;
cfg.compute.gpuRequired = false;
cfg.compute.resetGPU = false;

% -------------------------------------------------------------------------
% Robust LUT search
% -------------------------------------------------------------------------
cfg.lut.nGlobal = 120000;
cfg.lut.coarsePointCount = 450;
cfg.lut.globalFullRescoreTopK = 1800;

cfg.lut.nLocalSeeds = [8 5 3];
cfg.lut.nLocalPerSeed = [3000 2200 1400];
cfg.lut.localHalfWidth = [ ...
    20 20 15 1.0 1.0; ...
     8  8  6 0.5 0.5; ...
     3  3  2 0.2 0.2];

cfg.lut.localFullRescoreTopK = [1200 800 500];
cfg.lut.seedMinDistance = [0.80 0.65 0.50];
cfg.lut.peakSigmaQ_Ainv = 0.035;
cfg.lut.useGridIfSmall = false;

% Deterministic prior grid centered on the expected near-zero tilt geometry.
cfg.lut.priorGridXCount = 9;
cfg.lut.priorGridYCount = 9;
cfg.lut.priorGridZStep_mm = 5;
cfg.lut.priorMainTiltOffsets_deg = [-1 0 1];
cfg.lut.priorMinorTiltOffsets_deg = [-1 0 1];
cfg.lut.includeBroadAnchorGrid = true;
cfg.lut.broadAnchorXCount = 5;
cfg.lut.broadAnchorYCount = 5;
cfg.lut.broadAnchorZCount = 5;
cfg.lut.broadAnchorTiltStep_deg = 1;
cfg.lut.sobolSkipOffset = 0;

% -------------------------------------------------------------------------
% Ring assignment and exact refinement
% -------------------------------------------------------------------------
cfg.assign.maxDQ_Ainv = 0.075;
cfg.assign.minAssignedPixels = 80;
cfg.assign.maxRefinePixels = 5000;

cfg.opt.topKToRefine = 5;
cfg.opt.refineMinDistance = 0.18;
cfg.opt.reassignCycles = 3;
cfg.opt.reassignParamTol = 2e-4;
cfg.opt.robustScaleQ_Ainv = 0.010;
cfg.opt.maxIter = 160;
cfg.opt.maxFunEvals = 900;
cfg.opt.solver = 'auto';
cfg.opt.stepTolerance = 1e-6;
cfg.opt.functionTolerance = 1e-8;
cfg.opt.optimalityTolerance = 1e-6;

% -------------------------------------------------------------------------
% Automatic recovery
% -------------------------------------------------------------------------
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
cfg.search.expandFraction = [0.50 0.50 0.35 0.50 0.50];
cfg.search.expandParameters = [true true true true true];

% Keep automatic recovery physically reasonable for this experiment.
cfg.search.absoluteLimits = [ ...
   -100   200; ...       % beamX
   1450  1750; ...       % beamY
    180   350; ...       % z [mm]
     -8     8; ...       % thetaX [deg]
     -8     8];          % thetaY [deg]

% -------------------------------------------------------------------------
% Radial/profile output
% -------------------------------------------------------------------------
cfg.profile.qEdges_Ainv = 0.3:0.02:30.0;
cfg.profile.trimPct = 37;

cfg.tripol.qEdges = 0.3:0.025:30.0;
cfg.tripol.fixedPhiEdges = linspace(-pi,pi,361);
cfg.tripol.applyTotalCorr = true;
cfg.tripol.minPixPerBin = 1;
cfg.tripol.minPixPerRing = 1;
cfg.tripol.targetPixPerBin = 80;
cfg.tripol.binMean = 'trimmean';
cfg.tripol.trimPct = 37;

% -------------------------------------------------------------------------
% Output
% -------------------------------------------------------------------------
cfg.output.runTripol = false;
cfg.output.tripolEngine = 'fast';
cfg.output.makePlots = true;
cfg.output.saveMat = true;
cfg.output.matFile = 'lab6avg_tz400_pilatus3_2m_pose_fit.mat';

%% ========================================================================
%  RUN FIT
% =========================================================================
if exist(dataFile,'file') ~= 2
    error('Cannot find data file: %s',dataFile);
end
if exist('fit_lab6_pilatus3_2m_oldlut_robust.m','file') ~= 2
    error('Cannot find fit_lab6_pilatus3_2m_oldlut_robust.m on the MATLAB path.');
end

out = fit_lab6_pilatus3_2m_oldlut_robust(dataFile,cfg);

%% ========================================================================
%  PRINT RESULT
% =========================================================================
p = out.params;

fprintf('\n============================================================\n');
fprintf('FINAL PILATUS3 X 2M LaB6 POSE\n');
fprintf('============================================================\n');
fprintf('beamX_px   = %.6f\n',p(1));
fprintf('beamY_px   = %.6f\n',p(2));
fprintf('zBeam_mm   = %.6f\n',p(3));
fprintf('thetaX_deg = %.6f\n',p(4));
fprintf('thetaY_deg = %.6f\n',p(5));
fprintf('============================================================\n');
fprintf('Expected from direct image inspection: beam center near (49,1605) px,\n');
fprintf('z near 250-260 mm, and small tilt.\n');
fprintf('Total elapsed time: %.1f s\n',toc(tTotal));
