% Add dependencies.
addpath('..');
addpath(genpath('../common'));

%% Fiber-photometry data recorded with Doric DAQ.
inputDataFile = 'G:\My Drive\MSc\Research Data\Experiment 1.3 - FP and EMG Pilot\iso to wake\1L1R FP_2.csv';
% Columns corresponding to 465nm and 405nm.
signalColumn = 2;
referenceColumn = 4;
configuration = struct();
configuration.conditionEpochs = {'Pre', [0, 300], 'During', [300, 550], 'Post', [1010, 1310]};
configuration.bleachingEpochs = [-Inf, 300, 1010, Inf];
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