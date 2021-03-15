% Add dependencies.
addpath('..');
addpath(genpath('../common'));

% Fiber-photometry data recorded with Doric DAQ.
inputDataFile = '../data/Doric.csv';
% Columns corresponding to 465nm and 405nm.
signalColumn = 5;
referenceColumn = 3;
data = loadData(inputDataFile);
time = data(:, 1);
signal = data(:, signalColumn);
reference = data(:, referenceColumn);

configuration = struct();
configuration.conditionEpochs = {'Pre', [600, 700], 'During', [750, 850], 'Post', [900, 1000]};
configuration.baselineEpochs = [0, 300, 1010, Inf];
configuration.threshold = {2.91, @mad, @median};

% Call FPA with given configuration.
fpa = FPA(time, signal, reference, configuration);
cellfun(@warning, fpa.warnings);
fpa.plot();