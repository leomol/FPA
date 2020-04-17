% Add function dependencies.
addpath(genpath('../FPA'))

% Fiber-photometry data recorded with Doric DAQ.
fpFile = 'noclip2_2.csv';
% Columns corresponding to 465nm and 405nm.
fp465Column = 2;
fp405Column = 4;

% EMG data recorded with Axon.
emgFile = '19n28003V9.abf';
% Column corresponding to EMG (column 1 is time).
emgColumn = 4;

% Pre-processing configuration.
% EEG:
% Band-pass: 100-500
% abs --> rms
% Resample
% Smooth.
configuration = struct();
configuration.conditionEpochs = {'Awake', [700, 2000]};
configuration.bleachingEpochs = [700, 2000];
configuration.dffLowpassFrequency = 0.2;
configuration.peaksBandpassFrequency = [0.02, 0.2]; % 100/300
configuration.resamplingFrequency = 100;
configuration.f0Function = @movmean;
configuration.f0Window = 180;
configuration.f1Function = @movstd;
configuration.f1Window = 180;
configuration.thresholdingFunction = @mad;
configuration.thresholdFactor = 1;
configuration.plot = false;

% Process FP.
data = loadData(fpFile);
fpa = FPA(data(:, 1), data(:, fp465Column), data(:, fp405Column), configuration);
fpa.time = fpa.time(fpa.epochIds);
fpa.signal = fpa.signal(fpa.epochIds);
fpa.reference = fpa.reference(fpa.epochIds);
fpa.dff = fpa.dff(fpa.epochIds);
cellfun(@warning, fpa.warnings);
% Process EMG.
[data, units, names] = loadABF(emgFile);
emg = FPA(data(:, 1), data(:, emgColumn), [], configuration);
emg.time = emg.time(emg.epochIds);
emg.signal = emg.signal(emg.epochIds);
emg.reference = emg.reference(emg.epochIds);
emg.dff = emg.dff(emg.epochIds);
cellfun(@warning, emg.warnings);

% Compute envelope.
envelopeSamples = ceil(0.9 * configuration.resamplingFrequency);
[emg.up, emg.low] = envelope(emg.signal, envelopeSamples, 'rms');

% Plot raster.
figure('name', 'Raster');
hold('all');
plot(emg.time, emg.signal, 'DisplayName', 'EMG');
plot(emg.time, emg.up, 'LineWidth', 2, 'DisplayName', 'EMG envelope');
plot(fpa.time, fpa.signal, 'LineWidth', 2, 'DisplayName', 'FP');
title('FP and EMG');
xlabel('time (s)');
ylabel('z-score');
legend('show');

% Plot Cross-correlation.
figure('name', 'Cross-correlation');
xcorrSeconds = 60 * 10;
[xc, lags] = xcorr(fpa.signal, emg.up, xcorrSeconds * configuration.resamplingFrequency);
plot(lags / configuration.resamplingFrequency, xc);
title('Cross-correlation between FP and EMG');
xlabel('time (s)');
ylabel('xcorr');