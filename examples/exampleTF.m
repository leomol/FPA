% Add dependencies.
addpath('..');
addpath(genpath('../common'));

% Fiber-photometry data recorded with Doric DAQ.
inputDataFile = 'C:\Users\Molina\Documents\public\data\HALO\FibrePhotometry\Tamas\tf_m681_STS_d4_2.csv';
% Columns corresponding to 465nm and 405nm.
signalColumn = 2;
referenceColumn = 5;
data = loadData(inputDataFile);
time = data(:, 1);
signal = data(:, signalColumn);
reference = data(:, referenceColumn);

% General configuration to apply to all files.
configuration = struct();
configuration.resamplingFrequency = 100;
configuration.baselineLowpassFrequency = 0.1;
configuration.fitReference = true;
configuration.lowpassFrequency = 4;
configuration.peaksLowpassFrequency = 0.5;
configuration.threshold = {2.91, @mad, @median};
configuration.f0 = @mean;
configuration.f1 = @std;
configuration.conditionEpochs = {'Pre', [600, 700, 710, 720, 730, 740], 'During', [750, 850, 860, 870, 880, 890], 'Post', [900, 1000, 1010, 1020, 1030, 1040]};
configuration.baselineEpochs = [0, 300, 1010, Inf];

% Call FPA with given configuration.
fpa = FPA(time, signal, reference, configuration);
cellfun(@warning, fpa.warnings);
fpa.plot();
% fpa.export('exported');