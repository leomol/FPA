%% Add dependencies.
addpath('common');
addpath(genpath([getenv('USERPROFILE'), '/Documents/MATLAB/TDTMatlabSDK/TDTSDK']));

%% Example 1 - Analyze fiber-photometry data recorded with Doric DAQ.
inputDataFile = 'data/Doric photometry data.csv';
signalTitle = 'AIn-1 - Demodulated(Lock-In)';
referenceTitle = 'AIn-2 - Demodulated(Lock-In)';
configuration.resamplingFrequency = 20;
configuration.bleachingCorrectionEpochs = [-Inf, 600, 960, Inf];
configuration.zScoreEpochs = [-Inf, 600];
configuration.conditionEpochs = {'Pre', [100, 220], 'During 1', [650, 890], 'Post', [1480, 1600]};
configuration.triggeredWindow = 10;
configuration.f0Function = @movmean;
configuration.f0Window = 10;
configuration.f1Function = @movmean;
configuration.f1Window = 10;
configuration.peaksLowpassFrequency = 0.2;
configuration.thresholdingFunction = @mad;
configuration.thresholdFactor = 0.10;
[data, names] = loadData(inputDataFile);
s = ismember(names, signalTitle);
r = ismember(names, referenceTitle);
time = data(:, 1);
signal = data(:, s);
reference = data(:, r);
FPA(time, signal, reference, configuration);

%% Example 2 - Analyze fiber-photometry data recorded with Doric DAQ and behavioral data recorded with CleverSys.
% Fibre photometry recording file.
inputDataFile = 'data/Doric photometry data.csv';
% CleverSys event file in seconds and the name of the target sheet within.
inputEventFile = {'data/CleverSys event data.xlsx', 'Trial 1'};
% Names of columns corresponding to 465nm and 405nm.
signalTitle = 'AIn-1 - Demodulated(Lock-In)';
referenceTitle = 'AIn-2 - Demodulated(Lock-In)';
% Other settings.
configuration.resamplingFrequency = 20;
configuration.bleachingCorrectionEpochs = [-Inf, 600, 960, Inf];
configuration.artifactEpochs = [603, 620, 910, 915];
configuration.zScoreEpochs = [-Inf, 600];
configuration.triggeredWindow = 10;
configuration.f0Function = @movmean;
configuration.f0Window = 10;
configuration.f1Function = @movmean;
configuration.f1Window = 10;
configuration.peaksLowpassFrequency = 0.2;
configuration.thresholdingFunction = @mad;
configuration.thresholdFactor = 2;
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
fprintf(fid, '%.3f\n', time(results.peaksId));
fclose(fid);

%% Example 3 - Analyze fiber-photometry data recorded with TDT DAQ.
inputFolder = 'C:\Users\molina\Downloads\GP_PVN_13a-190531-122516';
signalTitle = 'Dv1A';
referenceTitle = 'Dv2A';
configuration.resamplingFrequency = 20;
configuration.bleachingCorrectionEpochs = [1, 1748];
configuration.zScoreEpochs = [-Inf, Inf];
configuration.conditionEpochs = {'Baseline', [1, 900], 'Test', [1102, 1702]};
configuration.triggeredWindow = 10;
configuration.f0Function = @movmean;
configuration.f0Window = 10;
configuration.f1Function = @movmean;
configuration.f1Window = 10;
configuration.peaksLowpassFrequency = 0.2;
configuration.thresholdingFunction = @mad;
configuration.thresholdFactor = 0.10;
data = loadTDT(inputFolder, {signalTitle, referenceTitle});
time = data(:, 1);
signal = data(:, 2);
reference = data(:, 3);
FPA(time, signal, reference, configuration);