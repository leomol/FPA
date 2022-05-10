% Add dependencies.
addpath('..');
addpath(genpath('../common'));

% Fiber-photometry data recorded with Doric DAQ.
inputDataFile = 'C:\Users\Molina\Documents\public\data\HALO\FibrePhotometry\Elizabeth\5_04_0007.csv';
inputEventFile = 'C:\Users\Molina\Documents\public\data\HALO\FibrePhotometry\Elizabeth\5_04.csv';
% Columns corresponding to 465nm and 405nm.
signalColumn = 2;
referenceColumn = 4;
data = loadData(inputDataFile);
time = data(:, 1);
signal = data(:, signalColumn);
reference = data(:, referenceColumn);

% General configuration to apply to all files.
configuration = struct();
configuration.resamplingFrequency = 100;
configuration.baselineLowpassFrequency = 0.1;
configuration.fitReference = false;
configuration.lowpassFrequency = 4;
configuration.peaksLowpassFrequency = 0.5;
configuration.threshold = {2.91, @mad, @median};
configuration.f0 = @mean;
configuration.f1 = @std;
configuration.conditionEpochs = loadBinaryStates(inputEventFile);
configuration.baselineEpochs = [0, 300, 640, 940];

% Call FPA with given configuration.
fpa = FPA(time, signal, reference, configuration);
cellfun(@warning, fpa.warnings);
%fpa.plotTrace();
fpa.export('exported');