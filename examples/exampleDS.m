% Add dependencies.
addpath('..');
addpath(genpath('../common'));

% Fiber-photometry data recorded with Doric DAQ.
inputDataFile = '../data/DoricDS.csv';
% Columns corresponding to 465nm and 405nm.
signalColumn = 2;
referenceColumn = 4;
configuration = struct();
configuration.conditionEpochs = {'Pre', [0, 300], 'During', [300, 600], 'Post', [600, 900]};
configuration.bleachingEpochs = [-Inf, Inf];
configuration.lowpassFrequency = 2;
configuration.peaksLowpassFrequency = 0.5;
configuration.thresholdingFunction = @mad;
configuration.thresholdFactor = 2.91;
% In the options below:
%   f: the calcium response after correcting for bleaching and motion artifacts.
%   f0: value calculated from neighbors around each value of f.
%   window: length of the window enclosing neighbors of f.
% Choose one according to your preference.
option = 1;
switch option
    case 1
        % CSM Opto's preference:
        %   z-score ==> (f - mean(f0)) / std(f0)
        %   This operation forces the signal to be centered at 0 and have 1 standard deviation tall.
        configuration.f0Function = @movmean;
        configuration.f0Window = Inf;
        configuration.f1Function = @movstd;
        configuration.f1Window = Inf;
    case 2
        % Doric_photom_analysis.m if plotting variable "normDat".
        configuration.f0 = 0;
        configuration.f1 = 1;
    case 3
        % Literature:
        %   df/f ==> (f - f0) / f0
        %   f0 is calculated sometimes with movmean other times with movmedian.
        configuration.f0Function = @movmedian;
        configuration.f0Window = Inf;
        configuration.f1Function = @movmad;
        configuration.f1Window = Inf;
end
data = loadData(inputDataFile);
time = data(:, 1);
signal = data(:, signalColumn);
reference = data(:, referenceColumn);
% Call FPA with given configuration.
results = FPA(time, signal, reference, configuration);
cellfun(@warning, results.warnings);