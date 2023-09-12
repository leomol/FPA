% Fiber-photometry data recorded with Doric DAQ.
inputDataFile = 'C:\Users\Molina\Documents\public\data\HALO\FibrePhotometry\Sabin\Males_Photo_Feb2023_0022.doric';

% Columns corresponding to 465nm and 405nm.
signalColumn = 2;
referenceColumn = 3;
data = loadH5(inputDataFile);
time = data(:, 1);
signal = data(:, signalColumn);
reference = data(:, referenceColumn);
baselineEpochs = [0, 600];

% General configuration.
configuration = struct();
configuration.resamplingFrequency = 100;
configuration.baselineLowpassFrequency = 0.1;
configuration.fitReference = false;
configuration.lowpassFrequency = 10.0;
configuration.peakSeparation = 0.5;
configuration.threshold = @(fpa, data) 2.0 * std(data);
configuration.peakDetectionMode = 'prominence';
configuration.peakWindow = 2.5;
% Use a portion of the data to compute modified z-score.
configuration.f0 = {@median, baselineEpochs};
configuration.f1 = {@mad, baselineEpochs};
configuration.conditionEpochs = {'Baseline 1', baselineEpochs, 'Stimulation', [600, 1800], 'Baseline 2', [1800, 2700]};
% Use a portion of the data to model and correct for photo-bleaching.
configuration.baselineEpochs = baselineEpochs;

% Call FPA with given configuration.
fpa = FPA(time, signal, reference, configuration);
cellfun(@warning, fpa.warnings);

% Plot.
fpa.plotTrace();

% Export trigger data for further analysis. Rerun selecting lower resamplingFrequency if needed.
% fpa.export('exported');