% Add dependencies.
addpath('..');
addpath(genpath('../common'));

%% Fiber-photometry data recorded with Doric DAQ.
inputDataFile = 'C:\Users\Molina\Documents\public\HALO\data\FibrePhotometry\Elizabeth\1L1R FP_2.csv';
% Columns corresponding to 465nm and 405nm.
signalColumn = 2;
referenceColumn = 4;
configuration = struct();
configuration.conditionEpochs = {'Pre', [0, 300], 'During', [300, 550], 'Post', [1010, 1310]};
configuration.bleachingEpochs = [-Inf, 300, 1010, Inf];
configuration.lowpassFrequency = 2;
configuration.peaksLowpassFrequency = 0.5;
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