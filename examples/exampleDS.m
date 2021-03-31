% Add dependencies.
addpath('..');
addpath(genpath('../common'));

% Fiber-photometry data recorded with Doric DAQ.
inputDataFile = '../data/Doric.csv';
% Columns corresponding to 465nm and 405nm.
signalColumn = 5;
referenceColumn = 3;
data = loadData(inputDataFile);
time = data(:, 1);
signal = data(:, signalColumn);
reference = data(:, referenceColumn);

configuration = struct();
configuration.conditionEpochs = {'Pre', [0, 300], 'During', [300, 600], 'Post', [600, 900]};
configuration.baselineEpochs = [-Inf, Inf];
configuration.lowpassFrequency = 2;
configuration.peaksLowpassFrequency = 0.5;

thresholdingOption = 3;
switch thresholdingOption
    case 1
        % 2.91 median absolute deviations from the median.
        configuration.threshold = {2.91, @mad, @median};
    case 2
        % 3.00 median absolute deviations from zero.
        configuration.threshold = {3.00, @mad, 0};
    case 3
        % 2.00 standtard deviations from the mean.
        configuration.threshold = {2.00, @std, @mean};
end

% In the options below:
%   f: the calcium response after baseline correction and motion artifacts.
%   f0: value calculated from neighbors around each value of f.
%   window: length of the window enclosing neighbors of f.
% Choose one according to your preference.
normalizationOption = 1;
switch normalizationOption
    case 1
        % "z-score" ==> (f - mean(f0)) / std(f0)
        configuration.f0 = @mean;
        configuration.f1 = @std;
    case 2
        % "altered z-score" ==> (f - median(f0)) / median(f0)
configuration.f0 = @median;
configuration.f1 = @mad;
    case 3
        % Doric_photom_analysis.m if plotting variable "normDat".
        configuration.f0 = 0;
        configuration.f1 = 1;
    case 4
        % df/f ==> (f - f0) / f0
        configuration.f0 = @mean;
        configuration.f1 = @mean;
end

% Call FPA with given configuration.
fpa = FPA(time, signal, reference, configuration);
cellfun(@warning, fpa.warnings);
fpa.plot();

% Export data.
[folder, basename] = fileparts(inputDataFile);
prefix = fullfile(folder, basename);
fpa.export(prefix);