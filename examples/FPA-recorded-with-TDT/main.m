% Add dependencies.
addpath(genpath('common'));

% Fiber-photometry data recorded with TDT DAQ.
inputFolder = 'data/RA_FST_320-200320-093026';
signalTitle = 'Dv1A';
referenceTitle = 'Dv2A';
configuration = struct();
configuration.conditionEpochs = {'Baseline', [300, 960], 'Stress 1', [960, 970], 'Stress 2', [1860, 1870], 'post-stress', [1870, 2410]}; %%add more epochs with: 'Name', [10000, 10001]
configuration.baselineEpochs = [-Inf, Inf];
configuration.bleachingEpochs = [30, 2500];
configuration.dffLowpassFrequency = 0.2;
configuration.peaksBandpassFrequency = [0.02, 0.2];
configuration.bleachingLowpassFrequency = 0.1;
configuration.resamplingFrequency = 20;
configuration.thresholdingFunction = @mad;
configuration.thresholdFactor = 1;

% Normally, movmean takes the mean over a moving window of a given size.
% If the window has infinite length, then this becomes the mean of everything.
% If f0 and f1 are both the mean of everything, then:
% df/f ==> (x - mean(f)) / mean(f)
% where
%   x is each point in your trace.
%   f is all points in your trace.
configuration.f0Function = @movmean;
configuration.f0Window = Inf;
configuration.f1Function = @movmean;
configuration.f1Window = Inf;

data = loadTDT(inputFolder, {signalTitle, referenceTitle});
time = data(:, 1);
signal = data(:, 2);
reference = data(:, 3);
results = FPA(time, signal, reference, configuration);