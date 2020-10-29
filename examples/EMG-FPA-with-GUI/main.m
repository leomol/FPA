% 2020-04-23. Leonardo Molina.
% 2020-05-07. Last modified.
% Draft #2 of a GUI for analyzing fiber-photometry and EMG data.
% See setupAnalysis and recalculate in GUI.m for details on analysis steps.

% Add FPA dependencies.
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '../..')));

% Global configuration.
configuration = struct();
configuration.resamplingFrequency = 100;
configuration.conditionEpochs = {'Awake', [800, 2000]};
configuration.xcorrSeconds = 5;

% Fiber-photometry configuration.
configuration.fpa = struct();
% Fiber-photometry data recorded with Doric DAQ.
configuration.fpa.file = 'C:/Users/molina/Documents/public/HALO/data/EMGFPA/noclip2_2.csv';
% Columns corresponding to 465nm and 405nm.
configuration.fpa.fp465Column = 2;
configuration.fpa.fp405Column = 4;
configuration.fpa.dffLowpassFrequency = 0.5;
configuration.fpa.f0Function = @movmean;
configuration.fpa.f0Window = 180;
configuration.fpa.f1Function = @movstd;
configuration.fpa.f1Window = 180;

% EMG configuration.
configuration.emg = struct();
% EMG data recorded with Axon.
configuration.emg.file = 'C:/Users/molina/Documents/public/HALO/data/EMGFPA/19n28003V9.abf';
% Column corresponding to EMG (column 1 is time).
configuration.emg.emgColumn = 4;
configuration.emg.bandpassFrequency = [100, 500];
configuration.emg.envelopeSize = 0.9;
configuration.emg.envelopeLowpassFrequency = 1.0;
configuration.emg.envelopeThreshold = 1.1;

% Run the analysis.
obj = GUI(configuration);

% Adjust settings in the console like this:
obj.envelopeSize = 0.9;
obj.envelopeLowpassFrequency = 1.0;
obj.envelopeThreshold = 1.1;