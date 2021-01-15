% Add dependencies.
addpath('..');
addpath(genpath('../common'));

%% Example 1 - Fiber-photometry data recorded with Doric DAQ.
inputDataFile = '../data/Doric.csv';
% Columns corresponding to 465nm and 405nm.
signalColumn = 5;
referenceColumn = 3;
configuration = struct();
configuration.thresholdingFunction = @mad;
configuration.thresholdFactor = 2.91;
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
inputDataFile = '../data/Inscopix.csv';
inputEventFile = '../data/InscopixTTL.csv';
configuration = struct();
configuration.thresholdingFunction = @mad;
configuration.thresholdFactor = 2.91;
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
results = FPA(time, signal, [], configuration);
cellfun(@warning, results.warnings);

%% Example 3 - Fiber-photometry data recorded with Doric DAQ and behavioral data recorded with CleverSys.
% Fibre photometry recording file.
inputDataFile = '../data/Doric.csv';
% CleverSys event file in seconds and the name of the target sheet within.
inputEventFile = {'../data/CleverSys.xlsx', 'Trial 1'};
signalColumn = 5;
referenceColumn = 3;
configuration = struct();
configuration.thresholdingFunction = @mad;
configuration.thresholdFactor = 2.91;
% Extract epochs from CleverSys output.
configuration.conditionEpochs = loadCleverSys(inputEventFile{:});
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
inputDataFile1 = '../data/Doric.csv';
inputDataFile2 = '../data/Doric2.csv';
% Columns corresponding to 465nm and 405nm.
signalColumn = 5;
referenceColumn = 3;
configuration = struct();
configuration.thresholdingFunction = @mad;
configuration.thresholdFactor = 2.91;
configuration.plot = false;
% Parse file with baseline data.
data = loadData(inputDataFile1);
time = data(:, 1);
signal = data(:, signalColumn);
reference = data(:, referenceColumn);
baseline = FPA(time, signal, reference, configuration);
% Parse file with test data.
data = loadData(inputDataFile2);
time = data(:, 1);
signal = data(:, signalColumn);
reference = data(:, referenceColumn);
configuration.f0 = baseline.f0;
configuration.f1 = baseline.f1;
configuration.plot = true;
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
inputFolder = '../data/GP_PVN_13a-190531-122516';
signalTitle = 'Dv1A';
referenceTitle = 'Dv2A';
configuration = struct();
configuration.thresholdingFunction = @mad;
configuration.thresholdFactor = 2.91;
data = loadTDT(inputFolder, {signalTitle, referenceTitle});
time = data(:, 1);
signal = data(:, 2);
reference = data(:, 3);
FPA(time, signal, reference, configuration);
cellfun(@warning, results.warnings);

%% Example 6 - Fiber-photometry data recorded with Doric DAQ. Two conditions are encoded as TTL values 1 or 0.
inputDataFile = '../data/Doric.csv';
% Columns corresponding to 465nm and 405nm.
signalColumn = 5;
referenceColumn = 3;
ttlColumn = 6;
homeCageEpoch = [0, 10];
configuration = struct();
configuration.thresholdingFunction = @mad;
configuration.thresholdFactor = 2.91;
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