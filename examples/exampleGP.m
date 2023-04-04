% Fiber-photometry data recorded with Doric DAQ.
inputDataFile = 'C:\Users\Molina\Documents\public\data\HALO\FibrePhotometry\Gavin\202211/GP_JZLPho_JZL_Nov20_M43.csv';

% Columns corresponding to 465nm and 405nm.
signalColumn = 2;
referenceColumn = 4;
data = loadData(inputDataFile);
time = data(:, 1);
signal = data(:, signalColumn);
reference = data(:, referenceColumn);

baseline = [1100, 1280, 2700, 2880];

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
configuration.f0 = {@median, baseline};
configuration.f1 = {@mad, baseline};
configuration.conditionEpochs = {'Pre footshock 1', [1200, 1500], 'Footshock 1', [1600, 1900], 'Post footshock 1', [2000, 2300], 'Pre footshock 2', [8500, 8800], 'Footshock 2', [8900, 9200], 'Post footshock 2', [9300, 9600]};
% Use a portion of the data to model and correct for photo-bleaching.
configuration.baselineEpochs = baseline;

% Call FPA with given configuration.
fpa = FPA(time, signal, reference, configuration);
cellfun(@warning, fpa.warnings);

% Plot.
fpa.plotTrace();

% Export trigger data for further analysis. Rerun selecting lower resamplingFrequency if needed.
% fpa.export('exported');