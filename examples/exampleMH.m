% Add dependencies.
addpath('..');
addpath(genpath('../common'));

%% Fiber-photometry data recorded with Doric DAQ.
inputDataFile = 'data/20201104_LH_K_M46.csv';
signalColumn = 2;
configuration = struct();
configuration.conditionEpochs = {'Data', [2 * 60, 4 * 60]};
configuration.bleachingEpochs = [2 * 60, 4 * 60];
configuration.lowpassFrequency = 2;
configuration.peaksLowpassFrequency = 0.5;
configuration.thresholdingFunction = @mad;
configuration.thresholdFactor = 2.91;
data = loadData(inputDataFile);
time = data(:, 1);
signal = data(:, signalColumn);
reference = []; % There is no reference data.
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