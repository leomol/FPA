% Large data files exported from LabChart require an SDK for loading into MATLAB (it is not the typical memory limitation).
% Using Doric's demodulated-FP
% Ignoring Powerlab's raw-FP data.
% Behavior is aligned with FP and EMG using the 20 TTL pulses occurring at 600s.
% Change recording setup to save camera frame triggers.
% 
% FP1/FP2/EMG: continuosly from the start.
% Camera trigger: Once at the start.
% TTL sync input: 10min after start, Square wave (20 pulses, for 9.5s)

sessionFolder = 'C:/Users/Molina/Documents/public/HALO/data/EMGFPA/20200701/';
getFirst = @(list) fullfile(list(1).folder, list(1).name);
first = @(expression) getFirst(dir(fullfile(sessionFolder, expression)));

behaviorDuration = 1800;
resamplingFrequency = 100;
xcorrSeconds = 5;

emg = struct();
emg.envelopeSize = 0.9;
emg.envelopeLowpassFrequency = 1.0;
emg.bandpassFrequency = [20, 400];

fpa = struct();
fpa.signalColumn = 2;
fpa.referenceColumn = 4;
fpa.cameraColumn = 9;
fpa.bleachingEpochs = [-Inf, Inf];
fpa.peaksBandpassFrequency = [0.02, 0.2];
fpa.dffLowpassFrequency = fpa.peaksBandpassFrequency(2);
fpa.bleachingLowpassFrequency = 0.1;
fpa.resamplingFrequency = resamplingFrequency;
fpa.thresholdingFunction = @mad;
fpa.thresholdFactor = 2;
fpa.f0Function = @movmean;
fpa.f0Window = 120;
fpa.f1Function = @movstd;
fpa.f1Window = 120;
fpa.file = first('*.csv');

cleversys = struct();
cleversys.trial = 1;
cleversys.file = first('*.xlsx');

powerlab.fp465Channel = 1;
powerlab.fp405Channel = 2;
powerlab.emgChannel = 3;
powerlab.cameraChannel = 6;
powerlab.file = first('*.adicht');

fprintf('Loading Doric data ...\n');
fpa.data = loadData(fpa.file);
fpa.time = fpa.data(:, 1);
fpa.signal = fpa.data(:, fpa.signalColumn);
fpa.reference = fpa.data(:, fpa.referenceColumn);
behaviorStart = fpa.time(find(fpa.data(:, fpa.cameraColumn), 1));

fprintf('Loading CleverSys data ...\n');
epochs = loadCleverSys(cleversys.file, cleversys.trial);
epochs(2:2:end) = cellfun(@(epoch) epoch + behaviorStart, epochs(2:2:end), 'UniformOutput', false);
behaviorEnd = max(cat(2, epochs{2:2:end}, behaviorStart + behaviorDuration));
preBaseline = [-Inf, behaviorStart];
postBaseline = [behaviorEnd, Inf];
epochs = ['pre-baseline', preBaseline, epochs, 'post-baseline', postBaseline];

fprintf('Loading LabChart data ...\n');
[time, data, units, names] = loadAdicht(powerlab.file);
time = time(:, end);
data = data(:, end);
cameraTime = time{powerlab.cameraChannel};
cameraData = data{powerlab.cameraChannel};
behaviorStart = cameraTime(find(cameraData > 0.5, 1));
emg.time = time{powerlab.emgChannel};
emg.signal = data{powerlab.emgChannel};

fprintf('Processing data ...\n');
[mn1, mx1] = bounds(fpa.time);
[mn2, mx2] = bounds(emg.time);
mn = max(mn1, mn2);
mx = min(mx1, mx2);
k1 = fpa.time >= mn & fpa.time <= mx;
fpa.time = fpa.time(k1);
fpa.signal = fpa.signal(k1);
fpa.reference = fpa.reference(k1);
k2 = emg.time >= mn & emg.time <= mx;
emg.time = emg.time(k2);
emg.signal = emg.signal(k2);

fpa.conditionEpochs = epochs;
fpa.plot = false;
[fpa.time, fpa.signal, fpa.reference] = alignTimestamps(fpa.time, 1 / resamplingFrequency, fpa.signal, fpa.reference);
fpa = FPA(fpa.time, fpa.signal, fpa.reference, fpa);
cellfun(@warning, fpa.warnings);

[emg.time, emg.signal] = alignTimestamps(emg.time, 1 / resamplingFrequency, emg.signal);
sourceFrequency = 1 / median(diff(emg.time));

% Option 1.
bandpassFilter = designfilt('bandpassiir', 'HalfPowerFrequency1', emg.bandpassFrequency(1), 'HalfPowerFrequency2', emg.bandpassFrequency(2), 'SampleRate', sourceFrequency, 'DesignMethod', 'butter', 'FilterOrder', 12);
emg.signal1 = filtfilt(bandpassFilter, emg.signal);

% Option 2.
bandpassFilter = designfilt('bandpassfir', 'CutoffFrequency1', emg.bandpassFrequency(1), 'CutoffFrequency2', emg.bandpassFrequency(2), 'SampleRate', sourceFrequency, 'DesignMethod', 'window', 'FilterOrder', 12);
emg.signal2 = filtfilt(bandpassFilter, emg.signal);

% Option 3.
fn = sourceFrequency / 2;
fl = emg.bandpassFrequency(1) / fn;
fh = emg.bandpassFrequency(2) / fn;
[b, a] = butter(4, [fl, fh], 'bandpass');
emg.signal3 = filter(b, a, emg.signal);

notchFilter = designfilt('bandstopiir', 'FilterOrder', 2, 'HalfPowerFrequency1', 59, 'HalfPowerFrequency2', 61, 'DesignMethod', 'butter', 'SampleRate', sourceFrequency);

emg.signal = filtfilt(notchFilter, emg.signal);
emg.signal = detrend(emg.signal);
if resamplingFrequency < sourceFrequency
    [p, q] = rat(resamplingFrequency / sourceFrequency);
    [emg.signal, emg.time] = resample(emg.signal, emg.time, resamplingFrequency, p, q);
else
    warning('[emg] Cannot resample to frequencies higher than the source frequency (%.2f Hz).', sourceFrequency);
end
emg.envelopeSamples = max(2, ceil(emg.envelopeSize * resamplingFrequency));
[high, low] = envelope(emg.signal, emg.envelopeSamples, 'rms');
lowpassFilter = designfilt('lowpassiir', 'HalfPowerFrequency', emg.envelopeLowpassFrequency, 'SampleRate', resamplingFrequency, 'DesignMethod', 'butter', 'FilterOrder', 12);
emg.high = filtfilt(lowpassFilter, high);
emg.low = filtfilt(lowpassFilter, low);
%%
fprintf('Plotting ...\n');
percentile = 0.99;
grow = 0.50;
cmap = lines();
color465 = [0.0000, 0.4470, 0.7410];
color405 = [0.8500, 0.3250, 0.0980];
colorDff = color465;
colorEmg = [1.0000, 0.2500, 0.2500];
colorEnvelope = [0.0000, 0.0000, 0.0000];
xlims = fpa.time([1, end]);
        
figureName = 'Raster plots';
figure('name', figureName);

subplot(3, 1, 1);
hold('all');
ylims = limits([fpa.signal; fpa.reference], percentile, grow);
plotEpochs(epochs, xlims, ylims, cmap, true);
plot(fpa.time, fpa.signal, 'Color', color465, 'DisplayName', '465');
plot(fpa.time, fpa.reference, 'Color', color405, 'DisplayName', '405');
legend('show');
ylim(ylims);

subplot(3, 1, 2);
hold('all');
ylims = limits(fpa.dff, percentile, grow);
plotEpochs(epochs, xlims, ylims, cmap, false);
plot(fpa.time, fpa.dff, 'Color', colorDff, 'DisplayName', 'df/f');
legend('show');
ylim(ylims);

subplot(3, 1, 3);
hold('all');
ylims = limits([emg.signal; emg.high], percentile, grow);
plotEpochs(epochs, xlims, ylims, cmap, false);
plot(emg.time, emg.signal, 'Color', colorEmg, 'DisplayName', 'EMG');
plot(emg.time, emg.high, 'Color', colorEnvelope, 'DisplayName', 'EMG envelope');
legend('show');
ylim(ylims);

linkaxes(findobj(gcf(), 'type', 'axes'), 'x');
xlim(xlims);

%%
figureName = 'Cross-correlation FP to EMG for different behaviors';
xcLags = round(xcorrSeconds * resamplingFrequency);
xcTics = (-xcLags:xcLags) / resamplingFrequency;
nEpochs = numel(epochs) / 2;
figure('name', figureName);
hold('all');
for e = 1:nEpochs
    epochName = epochs{2 * e - 1};
    ranges = epochs{2 * e};
    mask = time2id(fpa.time, ranges);
    xc = xcorr(fpa.dff(mask), emg.high(mask), xcLags, 'normalized');
    plot(xcTics, xc, 'DisplayName', epochName);
end
title(figureName);
legend('show');
xlabel('Lag (s)');

function plotEpochs(epochs, xlims, ylims, cmap, show)
    for e = 1:numel(epochs) / 2
        epochName = epochs{2 * e - 1};
        [faces, vertices] = patchEpochs(epochs{2 * e}, ylims(1), ylims(2));
        vertices(vertices == -inf) = xlims(1);
        vertices(vertices == +inf) = xlims(2);
        if show
            patch('Faces', faces, 'Vertices', vertices, 'FaceColor', cmap(e, :), 'EdgeColor', 'none', 'FaceAlpha', 0.50, 'DisplayName', sprintf('%s', epochName));
        else
            patch('Faces', faces, 'Vertices', vertices, 'FaceColor', cmap(e, :), 'EdgeColor', 'none', 'FaceAlpha', 0.50, 'HandleVisibility', 'off');
        end
    end
end

function [time2, varargout] = alignTimestamps(time, dt, varargin)
    % Shift points so that first time point is a multiple of dt.
    timestamp0 = ceil(time(1) / dt) * dt;
    time2 = colon(timestamp0, median(diff(time)), time(end))';
    varargout = cellfun(@(data) interp1(time(:), data(:), time2), varargin, 'UniformOutput', false);
end

function ylims = limits(x, percentile, grow)
    ylims = [prctile(x, 100 * (1 - percentile)), prctile(x, 100 * percentile)];
    delta = diff(ylims) * grow;
    ylims = [ylims(1) - delta, ylims(2) + delta];
end