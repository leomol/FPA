% Add dependencies.
addpath('..');
addpath(genpath('../common'));

% Fiber-photometry data recorded with TDT DAQ.
inputFolder = 'data/RA_FST_408-200813-112844';
signalTitle = 'Dv1A';
referenceTitle = 'Dv2A';
configuration = struct();
configuration.conditionEpochs = {'Baseline', [300, 900], 'Pickup 1', [900, 915], 'Intermission 1', [915, 1800], 'Pick up 2', [1800, 1815] 'Intermission 2', [1815, 2100] 'Swim', [2100, 3000] 'Post Stress', [3000, 3600]}; %%add more epochs with: 'Name', [10000, 10001]
configuration.bleachingEpochs = [30, 2500];
configuration.lowpassFrequency = 2;
configuration.peaksLowpassFrequency = 0.5;
configuration.thresholdingFunction = @mad;
configuration.thresholdFactor = 2.91;

data = loadTDT(inputFolder, {signalTitle, referenceTitle});
time = data(:, 1);
signal = data(:, 2);
reference = data(:, 3);
results = FPA(time, signal, reference, configuration);