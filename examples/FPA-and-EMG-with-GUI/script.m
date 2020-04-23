% Add function dependencies.
addpath(genpath('../FPA'));

% Global configuration.
resamplingFrequency = 100;
conditionEpochs = {'Awake', [800, 2000]};
xcorrSeconds = 5;

% Fiber-photometry analysis configuration.
fpa.configuration = struct();
% Fiber-photometry data recorded with Doric DAQ.
fpa.configuration.file = 'C:/Users/molina/Documents/public/MATLAB/EMGFPA/noclip2_2.csv';
% Columns corresponding to 465nm and 405nm.
fpa.configuration.fp465Column = 2;
fpa.configuration.fp405Column = 4;
% Signal processing settings.
fpa.configuration.bleachingEpochs = [700, 2000];
fpa.configuration.dffLowpassFrequency = 0.2;
fpa.configuration.f0Function = @movmean;
fpa.configuration.f0Window = 180;
fpa.configuration.f1Function = @movstd;
fpa.configuration.f1Window = 180;

% EMG configuration.
emg.configuration = struct();
% EMG data recorded with Axon.
emg.configuration.file = 'C:/Users/molina/Documents/public/MATLAB/EMGFPA/19n28003V9.abf';
% Column corresponding to EMG (column 1 is time).
emg.configuration.emgColumn = 4;
% Signal processing settings.
emg.configuration.bandpassFrequency = [100, 500];
emg.configuration.envelopeSize = 0.9;
emg.configuration.envelopeLowpassFrequency = 5;
emg.configuration.envelopeThreshold = 1;

% Complete configuration.
fpa.configuration.plot = false;
fpa.configuration.resamplingFrequency = resamplingFrequency;
fpa.configuration.conditionEpochs = conditionEpochs;
xcLags = round(xcorrSeconds * resamplingFrequency);

% Setup filters.
notchFilter = designfilt('bandstopiir', 'FilterOrder', 2, 'HalfPowerFrequency1', 59, 'HalfPowerFrequency2', 61, 'DesignMethod', 'butter', 'SampleRate', emg.sourceFrequency);
bandpassFilter = designfilt('bandpassfir', 'CutoffFrequency1', emg.configuration.bandpassFrequency(1), 'CutoffFrequency2', emg.configuration.bandpassFrequency(2), 'SampleRate', emg.sourceFrequency, 'DesignMethod', 'window', 'FilterOrder', 12);
lowpassFilter = designfilt('lowpassiir', 'HalfPowerFrequency', emg.configuration.envelopeLowpassFrequency, 'SampleRate', resamplingFrequency, 'DesignMethod', 'butter', 'FilterOrder', 12);

% Process FP.
% fpData = loadData(fpa.configuration.file);
fpa = FPA(fpData(:, 1), fpData(:, fpa.configuration.fp465Column), fpData(:, fpa.configuration.fp405Column), fpa.configuration);
fpa.time = fpa.time(fpa.epochIds);
fpa.signal = fpa.signal(fpa.epochIds);
fpa.reference = fpa.reference(fpa.epochIds);
fpa.dff = fpa.dff(fpa.epochIds);
% Print processing warnings.
cellfun(@warning, fpa.warnings);

% Process EMG.
emg.warnings = {};
% Band-pass (100 - 500 Hz), (rectify, resample, smooth ==> implemented with envelope function).
% [emgData, units, names] = loadABF(emg.configuration.file);
emg.time = emgData(:, 1);
emg.signal = emgData(:, emg.configuration.emgColumn);
emg.sourceFrequency = 1 / mean(diff(emg.time));
% Highpass filter.
emg.signal = filtfilt(bandpassFilter, emg.signal);
% Notch filter for 60 Hz.
emg.signal = filtfilt(notchFilter, emg.signal);
% Keep given time epochs.
ids = time2id(emg.time, [conditionEpochs{2:2:end}]);
emg.time = emg.time(ids);
emg.signal = emg.signal(ids);
% Detrend.
emg.signal = detrend(emg.signal);
% Resample to target frequency.
if resamplingFrequency < emg.sourceFrequency
    % Express frequency as a ratio p/q.
    [p, q] = rat(resamplingFrequency / emg.sourceFrequency);
    % Resample: interpolate every p/q/f, upsample by p, filter, downsample by q.
    [emg.signal, emg.time] = resample(emg.signal, emg.time, resamplingFrequency, p, q);
else
    emg.warnings{end + 1} = sprintf('[emg] Cannot resample to frequencies higher than the source frequency (%.2f Hz).', emg.sourceFrequency);
end

% Print processing warnings.
cellfun(@warning, emg.warnings);

% Envelope.
envelopeSamples = ceil(emg.configuration.envelopeSize * resamplingFrequency);
[high, low] = envelope(emg.signal, envelopeSamples, 'rms');
% Lowpass filter.
high = filtfilt(lowpassFilter, high);
low = filtfilt(lowpassFilter, low);
emg.envelope = [high, low];
% Threshold envelope.
emg.envelopeThreshold = mean(high) + std(high) * emg.configuration.envelopeThreshold;
emg.threshold = emg.envelopeThreshold * (high >= emg.envelopeThreshold);
% Cross-correlation.
xcTics = (-xcLags:xcLags) / resamplingFrequency;
[xc, lags] = xcorr(fpa.signal, emg.envelope(:, 1), xcLags);

% Initialize and format plot handles.
handles.isolatedFigure = figure('name', 'FPA and EMG');
handles.fpaAxes = subplot(2, 1, 1);
handles.emgAxes = subplot(2, 1, 2);
hold(handles.emgAxes, 'all');
linkaxes([handles.fpaAxes, handles.emgAxes], 'x');
handles.fpaPlot = plot(handles.fpaAxes, NaN(2, 1), NaN(2, 1), 'LineWidth', 1, 'DisplayName', 'Fiber-photometry');
handles.emgPlot = plot(handles.emgAxes, NaN(2, 1), NaN(2, 1), 'DisplayName', 'EMG');
handles.envelopePlot = plot(handles.emgAxes, NaN(2, 1), NaN(2, 1), 'LineWidth', 1, 'DisplayName', 'EMG envelope');
handles.thresholdPlot = plot(handles.emgAxes, NaN(2, 1), NaN(2, 1), 'LineWidth', 1, 'DisplayName', 'EMG threshold');
title(handles.fpaAxes, 'FP and EMG');
ylabel(handles.fpaAxes, 'z-score');
xlabel(handles.emgAxes, 'time (s)');
legend(handles.fpaAxes, 'show');
legend(handles.emgAxes, 'show');

% Update plots.
set(handles.fpaPlot, 'XData', fpa.time, 'YData', fpa.signal);
set(handles.emgPlot, 'XData', emg.time, 'YData', emg.signal);
set(handles.envelopePlot, 'XData', emg.time, 'YData', emg.envelope(:, 1));
set(handles.thresholdPlot, 'XData', emg.time, 'YData', emg.threshold);

% Cross-correlation plot.
handles.mixedFigure = figure('name', 'FPA vs EMG');
handles.xcAxes = axes();
handles.xcPlot = plot(handles.xcAxes, NaN(2, 1), NaN(2, 1), 'DisplayName', 'Cross-correlation');
title(handles.xcAxes, 'Cross-correlation between FP and EMG');
xlabel(handles.xcAxes, 'lag (s)');
ylabel(handles.xcAxes, 'xcorr');
set(handles.xcPlot, 'XData', xcTics, 'YData', xc);