%Robust high-angle PILATUS 100K LaB6 pose fit.
clear; clc; close all
tic

cfg = struct();

% -------------------------------------------------------------------------
% Detector / beam
% -------------------------------------------------------------------------
cfg.detector.pixel_um = 172;
cfg.beam.EkeV = 16.6;
cfg.beam.polMode = 'horizontal';
cfg.image.orientation = 'as_loaded';

% -------------------------------------------------------------------------
% Broad physical bounds
%
% At a large detector tilt, the direct-beam intercept can move far outside
% the active panel even when the detector center has moved only modestly.
% Therefore beamX must be allowed to be substantially negative.
% -------------------------------------------------------------------------
cfg.bounds.zBeam_mm = [70 250];
cfg.bounds.beamX_px = [-300 150];
cfg.bounds.beamY_px = [20 175];
cfg.bounds.thetaX_deg = [-15 15];
cfg.bounds.thetaY_deg = [0 70];

% This is a guide, not a hidden hard restriction.
cfg.prior.mainTiltAxis = 'y';
cfg.prior.tiltDeg = 30;
cfg.prior.tiltSign = +1;
cfg.prior.constrainBounds = false;
cfg.prior.tiltSigma_deg = 15;
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
% LUT scoring may use single precision; every promising candidate and the
% nonlinear fit are rescored/refined in exact double precision on the CPU.
% -------------------------------------------------------------------------
cfg.compute.device = 'gpu';             % 'auto', 'cpu', or 'gpu'
cfg.compute.autoBenchmark = true;
cfg.compute.gpuPrecision = 'single';
cfg.compute.cpuBlockCandidates = 512;
cfg.compute.gpuBlockCandidates = 4096;
cfg.compute.cpuParallel = 'auto';
cfg.compute.hitLookupDQ_Ainv = 1e-5;

% -------------------------------------------------------------------------
% Broad global LUT + three successively narrower local LUT levels
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
% Automatic recovery when a solution sticks to a bound or remains poor.
% For example, a fitted beamX equal to the lower bound triggers a new round
% with the beamX range expanded to the left.
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
% Reference / optional downstream output
% -------------------------------------------------------------------------
cfg.ref.cifFile = '1000057.cif';
cfg.output.runTripol = false;
cfg.output.tripolEngine = 'fast';
cfg.output.makePlots = true;
cfg.output.saveMat = true;
cfg.output.matFile = 'scan2_pilatus100k_oldlut_robust_pose_fit.mat';

out = fit_pilatus100k_lut_robust('scan2_0004.mat', cfg);
toc
