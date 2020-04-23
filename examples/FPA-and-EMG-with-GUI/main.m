% Draft of a GUI for analyzing fiber-photometry and EMG data.
% See analysis.m for details on analysis steps.
% See script.m for a non-GUI simpler version independent of main.m, analysis.m, and, Plots.m.

% Add function dependencies.
addpath(genpath('../FPA/../'));

% Global configuration.
shared.resamplingFrequency = 100;
shared.conditionEpochs = {'Awake', [800, 2000]};
shared.xcorrSeconds = 5;

% Fiber-photometry configuration.
fpa.configuration = struct();
% Fiber-photometry data recorded with Doric DAQ.
fpa.configuration.file = 'C:/Users/molina/Documents/public/HALO/data/EMGFPA/noclip2_2.csv';
% Columns corresponding to 465nm and 405nm.
fpa.configuration.fp465Column = 2;
fpa.configuration.fp405Column = 4;
fpa.configuration.bleachingEpochs = [700, 2000];
fpa.configuration.dffLowpassFrequency = 0.2;
fpa.configuration.f0Function = @movmean;
fpa.configuration.f0Window = 180;
fpa.configuration.f1Function = @movstd;
fpa.configuration.f1Window = 180;

% EMG configuration.
emg.configuration = struct();
% EMG data recorded with Axon.
emg.configuration.file = 'C:/Users/molina/Documents/public/HALO/data\EMGFPA/19n28003V9.abf';
% Column corresponding to EMG (column 1 is time).
emg.configuration.emgColumn = 4;
emg.configuration.bandpassFrequency = [100, 500];
emg.configuration.envelopeSize = 0.9;
emg.configuration.envelopeLowpassFrequency = 5;
emg.configuration.envelopeThreshold = 1;

% Run the analysis.
analysis(shared, fpa, emg);