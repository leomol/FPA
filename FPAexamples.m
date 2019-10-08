%% Add dependencies.
addpath(genpath('common'));

%% Example 1 - Fiber-photometry data recorded with Doric DAQ.
inputDataFile = 'data/Doric.csv';
% Names of columns corresponding to 465nm and 405nm.
signalTitle = 'AIn-1 - Demodulated(Lock-In)';
referenceTitle = 'AIn-2 - Demodulated(Lock-In)';
% Configuration (help FPA).
configuration = struct();
configuration.conditionEpochs = {'Pre', [100, 220], 'During', [650, 890], 'Post', [1480, 1600]};
configuration.bleachingEpochs = [-Inf, 600, 960, Inf];
configuration.dffLowpassFrequency = 0.2;
configuration.peaksBandpassFrequency = [0.02, 0.2];
configuration.bleachingLowpassFrequency = 0.1;
configuration.resamplingFrequency = 100;
configuration.f0Function = @movmean;
configuration.f0Window = 600;
configuration.f1Function = @movmean;
configuration.f1Window = 600;
configuration.thresholdingFunction = @mad;
configuration.thresholdFactor = 0.1;
% Load data (help loadData).
[data, names] = loadData(inputDataFile);
% Identify columns in data.
s = ismember(names, signalTitle);
r = ismember(names, referenceTitle);
time = data(:, 1);
signal = data(:, s);
reference = data(:, r);
% Call FPA with given configuration.
FPA(time, signal, reference, configuration);

%% Example 2 - Fiber-photometry data with stimuli recorded with Inscopix.
inputDataFile = 'data/Inscopix.csv';
inputEventFile = 'data/InscopixTTL.csv';
configuration = struct();
configuration.f0Function = @movmean;
configuration.f0Window = 600;
configuration.f1Function = @movmean;
configuration.f1Window = 600;
configuration.thresholdingFunction = @mad;
configuration.thresholdFactor = 0.1;
configuration.triggeredWindow = 1;
configuration.dffLowpassFrequency = 2;
configuration.peaksBandpassFrequency = [0.2, 2];
data = loadData(inputDataFile);
time = data(:, 1);
signal = data(:, 2);
% Extract epochs from TTL output.
eventDuration = 3;
baselineOffset = -4;
events = loadInscopixTTL(inputEventFile);
events = [events; events + eventDuration];
epochs = {'Baseline', events + baselineOffset, 'Stimulation', events};
configuration.conditionEpochs = epochs;
FPA(time, signal, [], configuration);


%% Example 3 - Fiber-photometry data recorded with Doric DAQ and behavioral data recorded with CleverSys.
% Fibre photometry recording file.
inputDataFile = 'data/Doric.csv';
% CleverSys event file in seconds and the name of the target sheet within.
inputEventFile = {'data/CleverSys.xlsx', 'Trial 1'};
signalTitle = 'AIn-1 - Demodulated(Lock-In)';
referenceTitle = 'AIn-2 - Demodulated(Lock-In)';
configuration = struct();
configuration.resamplingFrequency = 20;
configuration.f0Function = @movmean;
configuration.f0Window = 60;
configuration.f1Function = @movmean;
configuration.f1Window = 60;
configuration.baselineEpochs = [-Inf, 600];
configuration.artifactEpochs = [603, 620, 910, 915];
configuration.bleachingEpochs = [-Inf, 600, 960, Inf];
configuration.thresholdingFunction = @mad;
configuration.thresholdFactor = 0.1;
% Extract epochs from CleverSys output.
events = loadCleverSysEvents(inputEventFile{:});
eventNames = events.keys;
configuration.conditionEpochs = cellfun(@(eventName) {eventName, reshape([events(eventName).start, events(eventName).start + events(eventName).duration]', 1, 2 * numel(events(eventName).start))}, eventNames, 'UniformOutput', false);
configuration.conditionEpochs = cat(2, configuration.conditionEpochs{:});
[data, names] = loadData(inputDataFile);
s = ismember(names, signalTitle);
r = ismember(names, referenceTitle);
time = data(:, 1);
signal = data(:, s);
reference = data(:, r);
results = FPA(time, signal, reference, configuration);
% Save peak times to file.
[folder, basename] = fileparts(inputDataFile);
output = fullfile(folder, sprintf('%s peak-time.csv', basename));
fid = fopen(output, 'w');
fprintf(fid, 'Peak Time (s)\n');
fprintf(fid, '%.3f\n', results.time(results.peaksId));
fclose(fid);

%% Example 4 - Fiber-photometry data recorded with TDT DAQ.
inputFolder = 'data/GP_PVN_13a-190531-122516';
signalTitle = 'Dv1A';
referenceTitle = 'Dv2A';
configuration = struct();
configuration.conditionEpochs = {'Baseline', [1, 900], 'Test', [1102, 1702]};
configuration.baselineEpochs = [-Inf, Inf];
configuration.bleachingEpochs = [1, 1748];
configuration.dffLowpassFrequency = 0.2;
configuration.peaksBandpassFrequency = [0.02, 0.2];
configuration.bleachingLowpassFrequency = 0.1;
configuration.resamplingFrequency = 20;
configuration.f0Function = @movmean;
configuration.f0Window = 10;
configuration.f1Function = @movmean;
configuration.f1Window = 10;
configuration.thresholdingFunction = @mad;
configuration.thresholdFactor = 0.10;
data = loadTDT(inputFolder, {signalTitle, referenceTitle});
time = data(:, 1);
signal = data(:, 2);
reference = data(:, 3);
FPA(time, signal, reference, configuration);