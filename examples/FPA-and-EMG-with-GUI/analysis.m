% 2020-04-23. Leonardo Molina.
% 2020-04-23. Last modified.
function analysis(shared, fpa, emg)
    % Complete configuration.
    fpa.configuration.plot = false;
    fpa.configuration.shared.resamplingFrequency = shared.resamplingFrequency;
    fpa.configuration.conditionEpochs = shared.conditionEpochs;
    fpa.data = loadData(fpa.configuration.file);
    emg.data = loadABF(emg.configuration.file);
    
    [shared, fpa, emg] = analysis1(shared, fpa, emg);
    [shared, fpa, emg] = analysis2(shared, fpa, emg);
    plots = Plots();
    plots.setData('fpa', fpa.time, fpa.signal, 'emg', emg.time, emg.signal, 'envelope', emg.time, emg.envelope(:, 1), 'threshold', emg.time, emg.threshold, 'xcorr', shared.xcorrTics, shared.xcorr);
    figs = findall(0, 'Type', 'figure');
    for fig = figs
        axis(findall(fig, 'Type', 'axes'), 'tight');
    end
    
    % GUI configuration.
    shared.envelopeSizeRange = [0, 3];
    shared.envelopeLowpassFrequencyRange = [1, shared.resamplingFrequency / 2];
    shared.envelopeThresholdRange = [0.01, 3];
    
    control = uifigure('name', sprintf('%s - control', mfilename('Class')), 'MenuBar', 'none', 'NumberTitle', 'off', 'ToolBar', 'none');
    
    dy = 50;
    nControls = 3;
    layout = uigridlayout(control);
    layout.RowHeight = repmat({dy}, [1, nControls]);
    layout.ColumnWidth = {'1x', '2x'};
    
    ui = uilabel(layout, 'Text', 'Envelope size:');
    ui.Layout.Column = 1;
    ui.Layout.Row = 1;
    value = min(max(emg.configuration.envelopeSize, shared.envelopeSizeRange(1)), shared.envelopeSizeRange(2));
    ui = uislider(layout, 'Limits', shared.envelopeSizeRange, 'Value', value, 'ValueChangedFcn', @(handle, event) update(handle, 'envelopeSize'));
    ui.Layout.Column = 2;
    ui.Layout.Row = 1;
    
    ui = uilabel(layout, 'Text', 'Envelope low-pass frequency:');
    ui.Layout.Column = 1;
    ui.Layout.Row = 2;
    value = min(max(emg.configuration.envelopeLowpassFrequency, shared.envelopeLowpassFrequencyRange(1)), shared.envelopeLowpassFrequencyRange(2));
    ui = uislider(layout, 'Limits', shared.envelopeLowpassFrequencyRange, 'Value', value, 'ValueChangedFcn', @(handle, event) update(handle, 'envelopeLowpassFrequency'));
    ui.Layout.Column = 2;
    ui.Layout.Row = 2;
    
    ui = uilabel(layout, 'Text', 'Envelope threshold:');
    ui.Layout.Column = 1;
    ui.Layout.Row = 3;
    value = min(max(emg.configuration.envelopeThreshold, shared.envelopeThresholdRange(1)), shared.envelopeThresholdRange(2));
    ui = uislider(layout, 'Limits', shared.envelopeThresholdRange, 'Value', value, 'ValueChangedFcn', @(handle, event) update(handle, 'envelopeThreshold'));
    ui.Layout.Column = 2;
    ui.Layout.Row = 3;
    
    position = control.Position;
    control.Position = [position(1:3), dy * (nControls + 1)];
    
    function update(handle, target)
        emg.configuration.(target) = handle.Value;
        [shared, fpa, emg] = analysis2(shared, fpa, emg);
        plots.setData('fpa', fpa.time, fpa.signal, 'emg', emg.time, emg.signal, 'envelope', emg.time, emg.envelope(:, 1), 'threshold', emg.time, emg.threshold, 'xcorr', shared.xcorrTics, shared.xcorr);
    end
end

function [shared, fpa, emg] = analysis1(shared, fpa, emg)
    % Process FP.
    result = FPA(fpa.data(:, 1), fpa.data(:, fpa.configuration.fp465Column), fpa.data(:, fpa.configuration.fp405Column), fpa.configuration);
    fnames = fieldnames(result);
    nNames = numel(fnames);
    for i = 1:nNames
        fname = fnames{i};
        fpa.(fname) = result.(fname);
    end
    
    fpa.time = fpa.time(fpa.epochIds);
    fpa.signal = fpa.signal(fpa.epochIds);
    fpa.reference = fpa.reference(fpa.epochIds);
    fpa.dff = fpa.dff(fpa.epochIds);
    % Print processing warnings.
    cellfun(@warning, fpa.warnings);

    % Process EMG.
    emg.warnings = {};
    % Band-pass (100 - 500 Hz), (rectify, resample, smooth ==> implemented with envelope function).
    emg.time = emg.data(:, 1);
    emg.signal = emg.data(:, emg.configuration.emgColumn);
    emg.sourceFrequency = 1 / mean(diff(emg.time));
    % Highpass filter.
    bandpassFilter = designfilt('bandpassfir', 'CutoffFrequency1', emg.configuration.bandpassFrequency(1), 'CutoffFrequency2', emg.configuration.bandpassFrequency(2), 'SampleRate', emg.sourceFrequency, 'DesignMethod', 'window', 'FilterOrder', 12);
    emg.signal = filtfilt(bandpassFilter, emg.signal);
    % Notch filter for 60 Hz.
    notchFilter = designfilt('bandstopiir', 'FilterOrder', 2, 'HalfPowerFrequency1', 59, 'HalfPowerFrequency2', 61, 'DesignMethod', 'butter', 'SampleRate', emg.sourceFrequency);
    emg.signal = filtfilt(notchFilter, emg.signal);
    % Keep given time epochs.
    ids = time2id(emg.time, [shared.conditionEpochs{2:2:end}]);
    emg.time = emg.time(ids);
    emg.signal = emg.signal(ids);
    % Detrend.
    emg.signal = detrend(emg.signal);
    % Resample to target frequency.
    if shared.resamplingFrequency < emg.sourceFrequency
        % Express frequency as a ratio p/q.
        [p, q] = rat(shared.resamplingFrequency / emg.sourceFrequency);
        % Resample: interpolate every p/q/f, upsample by p, filter, downsample by q.
        [emg.signal, emg.time] = resample(emg.signal, emg.time, shared.resamplingFrequency, p, q);
    else
        emg.warnings{end + 1} = sprintf('[emg] Cannot resample to frequencies higher than the source frequency (%.2f Hz).', emg.sourceFrequency);
    end
    
    % Print processing warnings.
    cellfun(@warning, emg.warnings);
end

function [shared, fpa, emg] = analysis2(shared, fpa, emg)
    % Envelope.
    envelopeSamples = max(2, ceil(emg.configuration.envelopeSize * shared.resamplingFrequency));
    [high, low] = envelope(emg.signal, envelopeSamples, 'rms');
    % Lowpass filter.
    lowpassFilter = designfilt('lowpassiir', 'HalfPowerFrequency', emg.configuration.envelopeLowpassFrequency, 'SampleRate', shared.resamplingFrequency, 'DesignMethod', 'butter', 'FilterOrder', 12);
    high = filtfilt(lowpassFilter, high);
    low = filtfilt(lowpassFilter, low);
    emg.envelope = [high, low];
    % Threshold envelope.
    emg.envelopeThreshold = mean(high) + std(high) * emg.configuration.envelopeThreshold;
    emg.threshold = emg.envelopeThreshold * (high >= emg.envelopeThreshold);
    % Cross-correlation.
    xcLags = round(shared.xcorrSeconds * shared.resamplingFrequency);
    shared.xcorrTics = (-xcLags:xcLags) / shared.resamplingFrequency;
    shared.xcorr = xcorr(fpa.signal, emg.envelope(:, 1), xcLags);
end