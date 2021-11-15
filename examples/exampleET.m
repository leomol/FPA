% Add dependencies.
addpath('..');
addpath(genpath('../common'));

% Fiber-photometry data recorded with Doric DAQ.
inputDataFile = 'Thesis\Experiment 1.5 - FP cohort\novel female 12-May-21\5.02_0003.csv';
% Columns corresponding to 465nm and 405nm.
signalColumn = 4;
referenceColumn = 2;
data = loadData(inputDataFile);
time = data(:, 1);
signal = data(:, signalColumn);
reference = data(:, referenceColumn);

configuration = struct();
configuration.conditionEpochs = {'Pre', [600, 700, 710, 720, 730, 740], 'During', [750, 850, 860, 870, 880, 890], 'Post', [900, 1000, 1010, 1020, 1030, 1040]};
configuration.baselineEpochs = [0, 300, 1010, Inf];
configuration.threshold = {2.91, @mad, @median};

% Call FPA with given configuration.
fpa = FPA(time, signal, reference, configuration);
cellfun(@warning, fpa.warnings);
fpa.plotTrace();