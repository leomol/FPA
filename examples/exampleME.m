% !!

% Add dependencies.
addpath('..');
addpath(genpath('../common'));

% Fiber-photometry data recorded with Deisseroth's Multifiber.
inputDataFile = 'data/example.csv';
% Columns corresponding to 465nm and 405nm.
configuration = struct();
configuration.conditionEpochs = {'Data', [10, Inf]};
configuration.bleachingEpochs = [10, Inf];
configuration.lowpassFrequency = 2;
configuration.peaksLowpassFrequency = 0.5;
configuration.bleachingLowpassFrequency = 0.1;
configuration.thresholdingFunction = @mad;
configuration.thresholdFactor = 2.91;
configuration.f0Function = @movmedian;
configuration.f0Window = Inf;
configuration.f1Function = @movstd;
configuration.f1Window = Inf;
configuration.plot = {'trace'};
data = readtable(inputDataFile);
frame = data.Frame_470nm;
time = (frame - 1) * 0.05;
signal = data.MeanInt_470nm;
reference = data.MeanInt_410nm;
% Call FPA with given configuration.
results = FPA(time, signal, reference, configuration);
cellfun(@warning, results.warnings);