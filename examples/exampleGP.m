% Add dependencies.
addpath('..');
addpath(genpath('../common'));

% Fiber-photometry data recorded with TDT DAQ.
inputFolder = '../data/TDT';
signalTitle = 'Dv1A';
referenceTitle = 'Dv2A';
data = loadTDT(inputFolder, {signalTitle, referenceTitle});
time = data(:, 1);
signal = data(:, 2);
reference = data(:, 3);

% Call FPA with given configuration.
configuration = struct();
configuration.conditionEpochs = {'Baseline', [300, 900], 'Pickup 1', [900, 915], 'Intermission 1', [915, 1800], 'Pick up 2', [1800, 1815] 'Intermission 2', [1815, 2100] 'Swim', [2100, 3000] 'Post Stress', [3000, 3600]}; %%add more epochs with: 'Name', [10000, 10001]
configuration.baselineEpochs = [30, 2500];
results = FPA(time, signal, reference, configuration);