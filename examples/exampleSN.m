% Fiber-photometry data recorded with Doric DAQ.
inputDataFile = 'C:\Users\Molina\Documents\MATLAB\Males_Photo_Feb2023_0022.doric';

% Columns corresponding to 465nm and 405nm.
signalColumn = 1;
referenceColumn = 2;
data = loadH5(inputDataFile, {
        '/Traces/Console/Time(s)/Console_time(s)'
        '/Traces/Console/AIn-1 - Dem (AOut-1)/AIn-1 - Dem (AOut-1)'
        '/Traces/Console/AIn-2 - Dem (AOut-2)/AIn-2 - Dem (AOut-2)'
    });
time = data(:, 1);
signal = data(:, signalColumn);
reference = data(:, referenceColumn);

% General configuration.
configuration = struct();
configuration.resamplingFrequency = 100;
configuration.baselineLowpassFrequency = 0.1;
configuration.fitReference = false;
configuration.lowpassFrequency = 10.0;
configuration.peakSeparation = 0.5;
configuration.threshold = @(data) 2.0 * std(data);
configuration.peakDetectionMode = 'prominence';
configuration.peakWindow = 2.5; % Change to -2.5:0.1:2.5 for lower resolution.
% Use a portion of the data to compute modified z-score.
%configuration.f0 = {@median, [-Inf, Inf]};
%configuration.f1 = {@mad, [-Inf, Inf]};
%configuration.conditionEpochs = {'Pre footshock 1', [1200, 1500], 'Footshock 1', [1600, 1900], 'Post footshock 1', [2000, 2300], 'Pre footshock 2', [8500, 8800], 'Footshock 2', [8900, 9200], 'Post footshock 2', [9300, 9600]};
% Use a portion of the data to model and correct for photo-bleaching.
%configuration.baselineEpochs = [1250, Inf];

% Call FPA with given configuration.
fpa = FPA(time, signal, reference, configuration);
cellfun(@warning, fpa.warnings);

% Plot.
fpa.plotTrace();

% Export trigger data for further analysis. Rerun selecting lower resamplingFrequency if needed.
% fpa.export('exported');