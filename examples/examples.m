% Add dependencies.
addpath('..');
addpath(genpath('../common'));

%% Example - Fiber-photometry data recorded with Doric DAQ.
inputDataFile = '../data/Doric.csv';
% Columns corresponding to 465nm and 405nm.
signalColumn = 5;
referenceColumn = 3;
data = loadData(inputDataFile);
time = data(:, 1);
signal = data(:, signalColumn);
reference = data(:, referenceColumn);
% Call FPA with given configuration.
configuration = struct();
FPA(time, signal, reference, configuration);

%% Example - Same as above with updated configuration.
inputDataFile = '../data/Doric.csv';
% Columns corresponding to 465nm and 405nm.
signalColumn = 5;
referenceColumn = 3;
data = loadData(inputDataFile);
time = data(:, 1);
signal = data(:, signalColumn);
reference = data(:, referenceColumn);

configuration = struct();
configuration.conditionEpochs = {'Data', [-Inf, Inf]};
configuration.artifactEpochs = [];
configuration.resamplingFrequency = 20;
configuration.baselineEpochs = [-Inf, Inf];
configuration.baselineLowpassFrequency = 0.1;
configuration.airPLS = false;
configuration.lowpassFrequency = 2;
configuration.peaksLowpassFrequency = 0.5;
configuration.fitReference = true;
configuration.f0 = @median;
configuration.f1 = @mad;
configuration.thresholdingFunction = @mad;
configuration.thresholdFactor = 2.91;
configuration.triggeredWindow = 10;
configuration.plot = {'trace', 'power', 'stats', 'trigger'};
% Call FPA with given configuration.
fpa = FPA(time, signal, reference, configuration);
cellfun(@warning, fpa.warnings);

%% Example - Fiber-photometry data with stimuli recorded with Inscopix.
inputDataFile = '../data/Inscopix.csv';
inputEventFile = '../data/InscopixTTL.csv';
data = loadData(inputDataFile);
time = data(:, 1);
signal = data(:, 2);
% Extract epochs from TTL output.
eventDuration = 3;
baselineOffset = -4;
events = loadInscopixTTL(inputEventFile);
events = [events, events + eventDuration]';
epochs = {'Baseline', events + baselineOffset, 'Stimulation', events};
configuration = struct();
configuration.conditionEpochs = epochs;
fpa = FPA(time, signal, [], configuration);
cellfun(@warning, fpa.warnings);

%% Example - Fiber-photometry data recorded with Doric DAQ. Two conditions are encoded as TTL values 1 or 0.
inputDataFile = '../data/Doric.csv';
% Columns corresponding to 465nm and 405nm.
signalColumn = 5;
referenceColumn = 3;
ttlColumn = 6;
homeCageEpoch = [0, 10];
[lowEpochs, highEpochs] = loadTTLEpochs(inputDataFile, ttlColumn);
lowEpochs(1) = max(homeCageEpoch(2), lowEpochs(1));
configuration = struct();
configuration.conditionEpochs = {'No stimulation', lowEpochs, 'Stimulation', highEpochs};
data = loadData(inputDataFile);
time = data(:, 1);
signal = data(:, signalColumn);
reference = data(:, referenceColumn);
% Call FPA with given configuration.
fpa = FPA(time, signal, reference, configuration);
cellfun(@warning, fpa.warnings);

%% Example - Fiber-photometry data recorded with Doric DAQ and behavioral data recorded with CleverSys.
% Fibre photometry recording file.
inputDataFile = '../data/Doric.csv';
% CleverSys event file in seconds and the name of the target sheet within.
inputEventFile = {'../data/CleverSys.xlsx', 'Trial 1'};
signalColumn = 5;
referenceColumn = 3;
data = loadData(inputDataFile);
time = data(:, 1);
signal = data(:, signalColumn);
reference = data(:, referenceColumn);
configuration = struct();
configuration.thresholdingFunction = @mad;
configuration.thresholdFactor = 2.91;
% Extract epochs from CleverSys output.
configuration.conditionEpochs = loadCleverSys(inputEventFile{:});
fpa = FPA(time, signal, reference, configuration);
cellfun(@warning, fpa.warnings);
% Save peak times to file.
[folder, basename] = fileparts(inputDataFile);
output = fullfile(folder, sprintf('%s peak-time.csv', basename));
fid = fopen(output, 'w');
fprintf(fid, 'Peak Time (s)\n');
fprintf(fid, '%.4f\n', fpa.time(fpa.peakIds));
fclose(fid);

%% Example - Fiber-photometry data recorded with Doric DAQ - baseline from another file.
inputDataFile1 = '../data/Doric.csv';
inputDataFile2 = '../data/Doric.csv';
% Columns corresponding to 465nm and 405nm.
signalColumn = 5;
referenceColumn = 3;
configuration = struct();
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
fpa = FPA(time, signal, reference, configuration);
cellfun(@warning, fpa.warnings);

%% Example - Fiber-photometry data recorded with TDT DAQ.
inputFolder = '../data/TDT';
signalTitle = 'Dv1A';
referenceTitle = 'Dv2A';
data = loadTDT(inputFolder, {signalTitle, referenceTitle});
time = data(:, 1);
signal = data(:, 2);
reference = data(:, 3);
configuration = struct();
FPA(time, signal, reference, configuration);
cellfun(@warning, fpa.warnings);