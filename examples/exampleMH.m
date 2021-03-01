% Add dependencies.
addpath('..');
addpath(genpath('../common'));

% Fiber-photometry data recorded with Doric DAQ.
inputDataFile = '../data/Doric.csv';
signalColumn = 5;
data = loadData(inputDataFile);
time = data(:, 1);
signal = data(:, signalColumn);
reference = []; % There is no reference channel.

configuration = struct();
configuration.conditionEpochs = {'Data', [2 * 60, 4 * 60]};
configuration.baselineEpochs = [2 * 60, 4 * 60];

% Call FPA with given configuration.
fpa = FPA(time, signal, reference, configuration);
cellfun(@warning, fpa.warnings);