% Add dependencies.
addpath('..');
addpath(genpath('../common'));

%% Example 1 - Fiber-photometry data recorded with Doric DAQ.
inputDataFile = 'data/Doric.csv';
% Columns corresponding to 465nm and 405nm.
signalColumn = 5;
referenceColumn = 3;
configuration = struct();
configuration.conditionEpochs = {'Pre', [100, 340], 'During', [650, 890], 'Post', [1200, 1440]};
configuration.bleachingEpochs = [-Inf, 600, 960, Inf];
configuration.dffLowpassFrequency = 0.2;
configuration.peaksBandpassFrequency = [0.02, 0.2];
configuration.bleachingLowpassFrequency = 0.1;
configuration.resamplingFrequency = 100;
configuration.f0Function = @movmean;
configuration.f0Window = 600;
configuration.f1Function = @movstd;
configuration.f1Window = 600;
configuration.thresholdingFunction = @mad;
configuration.thresholdFactor = 0.1;
data = loadData(inputDataFile);
time = data(:, 1);
signal = data(:, signalColumn);
reference = data(:, referenceColumn);
% Call FPA with given configuration.
results = FPA(time, signal, reference, configuration);
% Save dff for statistical analysis.
[folder, basename] = fileparts(inputDataFile);
output = fullfile(folder, sprintf('%s dff.csv', basename));
fid = fopen(output, 'w');
fprintf(fid, 'Time (s), df/f, epoch\n');
fprintf(fid, '%.4f, %.4f, %d\n', [results.time(results.epochIds), results.dff(results.epochIds), results.epochGroups]');
fclose(fid);
cellfun(@warning, results.warnings);

%% Example 2 - Fiber-photometry data with stimuli recorded with Inscopix.
inputDataFile = 'data/Inscopix.csv';
inputEventFile = 'data/InscopixTTL.csv';
configuration = struct();
configuration.f0Function = @movmean;
configuration.f0Window = 600;
configuration.f1Function = @movstd;
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
events = [events, events + eventDuration]';
epochs = {'Baseline', events + baselineOffset, 'Stimulation', events};
configuration.conditionEpochs = epochs;
FPA(time, signal, [], configuration);
cellfun(@warning, results.warnings);

%% Example 3 - Fiber-photometry data recorded with Doric DAQ and behavioral data recorded with CleverSys.
% Fibre photometry recording file.
inputDataFile = 'data/Doric.csv';
% CleverSys event file in seconds and the name of the target sheet within.
inputEventFile = {'data/CleverSys.xlsx', 'Trial 1'};
signalColumn = 5;
referenceColumn = 3;
configuration = struct();
configuration.resamplingFrequency = 20;
configuration.f0Function = @movmean;
configuration.f0Window = 60;
configuration.f1Function = @movstd;
configuration.f1Window = 60;
configuration.dffEpochs = [-Inf, 600];
configuration.artifactEpochs = [603, 620, 910, 915];
configuration.bleachingEpochs = [-Inf, 600, 960, Inf];
configuration.thresholdingFunction = @mad;
configuration.thresholdFactor = 0.1;
% Extract epochs from CleverSys output.
events = loadCleverSysEvents(inputEventFile{:});
eventNames = events.keys;
configuration.conditionEpochs = cellfun(@(eventName) {eventName, reshape([events(eventName).start, events(eventName).start + events(eventName).duration]', 1, 2 * numel(events(eventName).start))}, eventNames, 'UniformOutput', false);
configuration.conditionEpochs = cat(2, configuration.conditionEpochs{:});
data = loadData(inputDataFile);
time = data(:, 1);
signal = data(:, signalColumn);
reference = data(:, referenceColumn);
results = FPA(time, signal, reference, configuration);
% Save peak times to file.
[folder, basename] = fileparts(inputDataFile);
output = fullfile(folder, sprintf('%s peak-time.csv', basename));
fid = fopen(output, 'w');
fprintf(fid, 'Peak Time (s)\n');
fprintf(fid, '%.4f\n', results.time(results.peaksId));
fclose(fid);
% Save dff for statistical analysis.
output = fullfile(folder, sprintf('%s dff.csv', basename));
fid = fopen(output, 'w');
fprintf(fid, 'Time (s), df/f, epoch\n');
fprintf(fid, '%.4f, %.4f, %d\n', [results.time(results.epochIds), results.dff(results.epochIds), results.epochGroups]');
fclose(fid);
cellfun(@warning, results.warnings);

%% Example 4 - Fiber-photometry data recorded with Doric DAQ - baseline from another file.
inputDataFile1 = 'data/Doric.csv';
inputDataFile2 = 'data/Doric2.csv';
% Columns corresponding to 465nm and 405nm.
signalColumn = 5;
referenceColumn = 3;
configuration = struct();
configuration.dffLowpassFrequency = 0.2;
configuration.peaksBandpassFrequency = [0.02, 0.2];
configuration.bleachingLowpassFrequency = 0.1;
configuration.resamplingFrequency = 100;
configuration.f0Function = @movmean;
configuration.f0Window = Inf;
configuration.f1Function = @movstd;
configuration.f1Window = Inf;
configuration.thresholdingFunction = @mad;
configuration.thresholdFactor = 0.1;
% Parse file with baseline data.
data = loadData(inputDataFile1);
time = data(:, 1);
signal = data(:, signalColumn);
reference = data(:, referenceColumn);
baseline = FPA(time, signal, reference, configuration);
close(baseline.figures);
% Parse file with test data.
data = loadData(inputDataFile2);
time = data(:, 1);
signal = data(:, signalColumn);
reference = data(:, referenceColumn);
configuration.f0 = baseline.f0;
configuration.f1 = baseline.f1;
results = FPA(time, signal, reference, configuration);
% Save dff for statistical analysis.
[folder, basename] = fileparts(inputDataFile2);
output = fullfile(folder, sprintf('%s dff.csv', basename));
fid = fopen(output, 'w');
fprintf(fid, 'Time (s), df/f, epoch\n');
fprintf(fid, '%.4f, %.4f, %d\n', [results.time(results.epochIds), results.dff(results.epochIds), results.epochGroups]');
fclose(fid);
cellfun(@warning, results.warnings);

%% Example 5 - Fiber-photometry data recorded with TDT DAQ.
inputFolder = 'data/GP_PVN_13a-190531-122516';
signalTitle = 'Dv1A';
referenceTitle = 'Dv2A';
configuration = struct();
configuration.conditionEpochs = {'Baseline', [1, 900], 'Test', [1102, 1702]};
configuration.dffEpochs = [-Inf, Inf];
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
cellfun(@warning, results.warnings);

%% Example 6 - Fiber-photometry data recorded with Doric DAQ. Two conditions are encoded as TTL values 1 or 0.
inputDataFile = 'data/Doric.csv';
% Columns corresponding to 465nm and 405nm.
signalColumn = 5;
referenceColumn = 3;
ttlColumn = 6;
homeCageEpoch = [0, 10];
configuration = struct();
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
[lowEpochs, highEpochs] = loadTTLEpochs(inputDataFile, ttlColumn);
lowEpochs(1) = max(homeCageEpoch(2), lowEpochs(1));
configuration.conditionEpochs = {'No stimulation', lowEpochs, 'Stimulation', highEpochs};
data = loadData(inputDataFile);
time = data(:, 1);
signal = data(:, signalColumn);
reference = data(:, referenceColumn);
% Call FPA with given configuration.
results = FPA(time, signal, reference, configuration);
% Save dff for statistical analysis.
[folder, basename] = fileparts(inputDataFile);
output = fullfile(folder, sprintf('%s dff.csv', basename));
fid = fopen(output, 'w');
fprintf(fid, 'Time (s), df/f, epoch\n');
fprintf(fid, '%.4f, %.4f, %d\n', [results.time(results.epochIds), results.dff(results.epochIds), results.epochGroups]');
fclose(fid);
cellfun(@warning, results.warnings);

%% Example 7.
inputDataFile = 'C:\Users\molina\Documents\public\HALO\data\FibrePhotometry\Tamas\negative-peak-test.xlsx';
signalColumn = 10;
configuration = struct();
configuration.conditionEpochs = {'Data', [10, 170]};
configuration.bleachingEpochs = [10, 170];
configuration.dffEpochs = [10, 170];
configuration.dffLowpassFrequency = 2.00;
configuration.peaksBandpassFrequency = [0.02, 2.00];
configuration.bleachingLowpassFrequency = 0.1;
configuration.resamplingFrequency = 10;
configuration.f0Function = @movmean;
configuration.f0Window = 600;
configuration.f1Function = @movstd;
configuration.f1Window = 600;
configuration.movingWindow = false;
configuration.thresholdingFunction = @mad;
configuration.thresholdFactor = 1.5;
data = loadData(inputDataFile, 6);
time = data(:, 1);
signal = data(:, signalColumn);
reference = [];
% Call FPA with given configuration.
results = FPA(time, signal, reference, configuration);
cellfun(@warning, results.warnings);