% 2023-02-09. LM.
% 2023-02-09. Last modified.

filenames = {'data/M294_Exercise_Test_0001.doric', ...
             'data/M290_Exercise_Test_0001.doric', ...
            };
epochs = {...
           {'Baseline1', [0, 300], 'Pickup', [300, 480], 'Exercise', [480 4080], 'Baseline2', [4080, Inf]}, ...
           {'Baseline1', [0, 300], 'Pickup', [300, 480], 'Exercise', [480 4080], 'Baseline2', [4080, Inf]}, ...
         };

nFiles = numel(filenames);
for i = 1:nFiles
    filename = filenames{i};
    fprintf('Loading "%s"\n', filename);
    epoch = epochs{i};
    analysis(filename, epoch);
end
disp('Done!');

function analysis(filename, epoch)
    % Columns corresponding to 465nm and 405nm.
    signalColumn = 2;
    referenceColumn = 3;
    data = loadH5(filename);
    time = data(:, 1);
    signal = data(:, signalColumn);
    reference = data(:, referenceColumn);
    
    % General configuration.
    configuration = struct();
    configuration.resamplingFrequency = 100;
    configuration.baselineLowpassFrequency = 0.1;
    configuration.baselineFitType = 'exp1'; % exp1 | poly1
    configuration.lowpassFrequency = 10.0;
    configuration.peakSeparation = 0.5;
    configuration.threshold = @(fpa, data) 2.0 * std(data);
    configuration.peakDetectionMode = 'prominence';
    configuration.peakWindow = 2.5; % Change to -2.5:0.1:2.5 for lower resolution.
    % Use a portion of the data to compute modified z-score.
    range = [epoch{[2, 8]}];
    configuration.f0 = {@median, range};
    configuration.f1 = {@mad, range};
    configuration.conditionEpochs = epoch;
    % Use a portion of the data to model and correct for photo-bleaching.
    configuration.baselineEpochs = range;
    
    % Call FPA with given configuration.
    fpa = FPA(time, signal, reference, configuration);
    cellfun(@warning, fpa.warnings);
    
    fpa.export(filename);
    fpa.plotTrace();
    %fpa.plotStatistics();
    %fpa.plotTrigger();
    %fpa.plotTriggerAverage();
    %fpa.plotAUC();
end