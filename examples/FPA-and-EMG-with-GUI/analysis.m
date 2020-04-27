% 2020-04-23. Leonardo Molina.
% 2020-04-27. Last modified.
function analysis(shared, fpa, emg)
    % Complete configuration.
    fpa.configuration.plot = false;
    fpa.configuration.resamplingFrequency = shared.resamplingFrequency;
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
    shared.envelopeSizeRange = [0, 10];
    shared.envelopeLowpassFrequencyRange = [1e-6, shared.resamplingFrequency / 20];
    shared.envelopeThresholdRange = [1e-3, prctile(emg.signal, 99)];
    
    control = uifigure('name', sprintf('%s - control', mfilename('Class')), 'MenuBar', 'none', 'NumberTitle', 'off', 'ToolBar', 'none', 'CloseRequestFcn', @(handle, event)close());
    
    dy = 50;
    nControls = 3;
    layout = uigridlayout(control);
    layout.RowHeight = repmat({dy}, [1, nControls]);
    layout.ColumnWidth = {'2x', '3x'};
    
    target = 'envelopeSize';
    text = 'Envelope size';
    value = min(max(emg.configuration.envelopeSize, shared.envelopeSizeRange(1)), shared.envelopeSizeRange(2));
    label = uilabel(layout);
    label.Layout.Column = 1;
    label.Layout.Row = 1;
    slider = uislider(layout, 'Limits', shared.envelopeSizeRange, 'Value', value, 'ValueChangedFcn', @(slider, event) update(slider, label, text, target));
    slider.Layout.Column = 2;
    slider.Layout.Row = 1;
    update(slider, label, text, target);
    
    target = 'envelopeLowpassFrequency';
    text = 'Envelope low-pass frequency';
    value = min(max(emg.configuration.envelopeLowpassFrequency, shared.envelopeLowpassFrequencyRange(1)), shared.envelopeLowpassFrequencyRange(2));
    label = uilabel(layout);
    label.Layout.Column = 1;
    label.Layout.Row = 2;
    slider = uislider(layout, 'Limits', shared.envelopeLowpassFrequencyRange, 'Value', value, 'ValueChangedFcn', @(slider, event) update(slider, label, text, target));
    slider.Layout.Column = 2;
    slider.Layout.Row = 2;
    update(slider, label, text, target);
    
    target = 'envelopeThreshold';
    text = 'Envelope threshold';
    value = min(max(emg.configuration.envelopeThreshold, shared.envelopeThresholdRange(1)), shared.envelopeThresholdRange(2));
    label = uilabel(layout);
    label.Layout.Column = 1;
    label.Layout.Row = 3;
    slider = uislider(layout, 'Limits', shared.envelopeThresholdRange, 'Value', value, 'ValueChangedFcn', @(slider, event) update(slider, label, text, target));
    slider.Layout.Column = 2;
    slider.Layout.Row = 3;
    update(slider, label, text, target);
    
    position = control.Position;
    control.Position = [position(1:3), dy * (nControls + 1)];
    
    function update(slider, label, baseText, target)
        emg.configuration.(target) = slider.Value;
        label.Text = sprintf('%s (%.2f):', baseText, slider.Value);
        [shared, fpa, emg] = analysis2(shared, fpa, emg);
        plots.setData('fpa', fpa.time, fpa.signal, 'emg', emg.time, emg.signal, 'envelope', emg.time, emg.envelope(:, 1), 'threshold', emg.time, emg.threshold, 'xcorr', shared.xcorrTics, shared.xcorr);
    end
    
    function close()
        plots.close();
        delete(control);
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
    emg.envelopeThreshold = emg.configuration.envelopeThreshold;
    emg.threshold = prctile(emg.signal, 99) * (high >= emg.envelopeThreshold);
    % Cross-correlation.
    xcLags = round(shared.xcorrSeconds * shared.resamplingFrequency);
    shared.xcorrTics = (-xcLags:xcLags) / shared.resamplingFrequency;
    shared.xcorr = xcorr(fpa.signal, emg.envelope(:, 1), xcLags, 'normalized');
end