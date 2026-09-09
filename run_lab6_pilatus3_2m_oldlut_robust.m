%RUN_LAB6_PILATUS3_2M_OLDLUT_ROBUST
% One-command LaB6 pose fit for a PILATUS3 X 2M detector.
%
% Main fitted parameters:
%   p = [beamX_px beamY_px zBeam_mm thetaX_deg thetaY_deg]
%
% Required:
%   fit_lab6_pilatus3_2m_oldlut_robust.m
%   cif2powder_CM.m
%   1000057.cif   (LaB6)
%
% The detector raster is expected to be 1679 x 1475 [rows x cols], or the
% transposed orientation, with 172 um pixels.  Nominal 17-pixel vertical and
% 7-pixel horizontal inter-module gaps are masked explicitly.

clear; clc; close all
tic

% -------------------------------------------------------------------------
% EDIT THESE FIRST
% -------------------------------------------------------------------------
dataFile = 'lab6avg_tz400.mat';    % MAT file containing the 2-D detector image
EkeV = 59.04;                       % X-ray energy used for this LaB6 image

cfg = struct();

% -------------------------------------------------------------------------
% Detector / beam
% -------------------------------------------------------------------------
cfg.detector.name = 'PILATUS3_X_2M';
cfg.detector.pixel_um = 172;
cfg.detector.expectedSize_rc = [1679 1475];
cfg.detector.strictSize = false;   % true if you want a hard size check

cfg.beam.EkeV = EkeV;
cfg.beam.polMode = 'horizontal';
cfg.image.orientation = 'as_loaded';

% -------------------------------------------------------------------------
% Broad physical bounds
%
% Leave beamX/beamY empty for detector-size-derived broad bounds.  If you
% know the approximate direct-beam intercept, replacing these with tighter
% bounds is the single best way to accelerate the fit.
% -------------------------------------------------------------------------
% Direct beam is near / slightly outside the LOWER-LEFT corner
cfg.bounds.beamX_px = [-150  250];
cfg.bounds.beamY_px = [1450 1900];

 
% Detector expected nearly normal to beam
cfg.bounds.thetaX_deg = [-2 2];
cfg.bounds.thetaY_deg = [-2 2];

cfg.bounds.zBeam_mm = [150 400];




beamX_px     = [30 70];       % expected ~49.5
beamY_px     = [1580 1635];   % expected ~1607
zBeam_mm     = [245 275];     % expected ~258-259 mm
thetaX_deg   = [-1.5 1.5];    % expected ~0 deg
thetaY_deg   = [-1.5 1.5];    % expected ~+0.2 deg
EkeV         = 59.04;

% No 14-degree prior for this geometry
cfg.prior.tiltDeg = 0;
cfg.prior.mainTiltAxis = 'y';   % effectively irrelevant when tiltDeg = 0
cfg.prior.tiltSign = +1;



% Guide only, not a hidden hard restriction.
cfg.prior.mainTiltAxis = 'y';
cfg.prior.tiltDeg = 0;
cfg.prior.tiltSign = +1;
cfg.prior.constrainBounds = false;
cfg.prior.tiltSigma_deg = 0;
cfg.prior.tiltPenaltyWeight = 0.003;
cfg.prior.penaltyType = 'soft';
cfg.prior.maxPenalty = 0.025;

% -------------------------------------------------------------------------
% Signal processing
% -------------------------------------------------------------------------
cfg.prep.backgroundWindowPx = 35;
cfg.prep.thresholdSigma = 2.5;
cfg.prep.thresholdPercentile = 98.5;
cfg.prep.maxScorePixels = 6000;
cfg.prep.autoRelaxThreshold = true;
cfg.prep.minScorePixels = 1200;
cfg.prep.relaxSigmaFloor = 2.1;
cfg.prep.relaxPercentileFloor = 98.0;

% -------------------------------------------------------------------------
% CPU / GPU
% LUT scoring may use single precision.  Promising candidates and the final
% nonlinear refinement are rescored/refined in exact double precision.
% -------------------------------------------------------------------------
cfg.compute.device = 'auto';       % 'auto', 'cpu', or 'gpu'
cfg.compute.autoBenchmark = true;
cfg.compute.gpuPrecision = 'single';
cfg.compute.cpuBlockCandidates = 512;
cfg.compute.gpuBlockCandidates = 4096;
cfg.compute.cpuParallel = 'auto';
cfg.compute.hitLookupDQ_Ainv = 1e-5;

% -------------------------------------------------------------------------
% Same robust broad LUT + three local LUT levels as the latest 100K version
% -------------------------------------------------------------------------
cfg.lut.nGlobal = 120000;
cfg.lut.coarsePointCount = 450;
cfg.lut.globalFullRescoreTopK = 1800;
cfg.lut.nLocalSeeds = [8 5 3];
cfg.lut.nLocalPerSeed = [3000 2200 1400];
cfg.lut.localHalfWidth = [ ...
    45 32 28 12 14; ...
    18 14 12  5  6; ...
     6  6  4  2  2];
cfg.lut.localFullRescoreTopK = [1200 800 500];
cfg.lut.seedMinDistance = [0.80 0.65 0.50];
cfg.lut.peakSigmaQ_Ainv = 0.035;
cfg.lut.useGridIfSmall = false;

% -------------------------------------------------------------------------
% Automatic recovery if a solution sticks to a bound or remains poor
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
cfg.search.expandFraction = [1.0 0.60 0.60 0.50 0.50];

% -------------------------------------------------------------------------
% Exact final refinement
% -------------------------------------------------------------------------
cfg.opt.topKToRefine = 5;
cfg.opt.reassignCycles = 3;
cfg.opt.solver = 'auto';

% -------------------------------------------------------------------------
% LaB6 reference / output
% -------------------------------------------------------------------------
cfg.ref.cifFile = '1000057.cif';
cfg.output.runTripol = false;
cfg.output.tripolEngine = 'fast';
cfg.output.makePlots = true;
cfg.output.saveMat = true;
cfg.output.matFile = 'lab6_pilatus3_2m_oldlut_robust_pose_fit.mat';

out = fit_lab6_pilatus3_2m_oldlut_robust(dataFile, cfg);

toc
