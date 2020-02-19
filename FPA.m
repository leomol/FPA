% FPA(time, signal, reference, configuration)
% 
% Correct signal from bleaching and artifacts, normalize and detect peaks of
% spontaneous activity based on the given parameters.
% Signal and reference are column vectors.
% 
% Overall analysis steps:
% 	-Resample data.
% 	-Fit an exponential decay in a portion of the data representative of bleaching (low-pass).
% 	-Normalize signal to reference.
% 	-Compute df/f in a (moving) window.
% 	-Find peaks of spontaneous activity (band-pass).
% 	-Plot corrected signal and peaks; highlight epochs.
% 	-Plot triggered averages.
% 	-Plot power spectrum.
% 	-Plot stats.
% Note that some plots are of df/f while others are of filtered df/f.
% 
% configuration is a struct with the following fields (defaults are used for missing fields):
%     conditionEpochs - Epochs for different conditions: {'epoch1', [start1, end1, start2, end2, ...], 'epoch2', ...}
%     dffEpochs - Time epochs (s) to include for df/f normalization.
%     bleachingEpochs - Time epochs (s) to include for bleaching correction.
%     artifactEpochs - Time epochs (s) to remove.
%     resamplingFrequency - Resampling frequency (Hz).
%     dffLowpassFrequency - Lowpass frequency to filter df/f.
%     peaksBandpassFrequency - Low/High frequencies to detect peaks.
%     bleachingLowpassFrequency - Lowpass frequency to detect bleaching decay.
%     f0Function - One of @movmean, @movmedian, @movmin.
%     f0Window - Length of the moving window to calculate f0.
%     f1Function - One of @movmean, @movmedian, @movmin, @movstd.
%     f1Window - Length of the moving window to calculate f1.
%     thresholdingFunction - One of @mad, @std.
%     thresholdFactor - Thresholding cut-off.
%     triggeredWindow - Length of time to capture around each peak of spontaneous activity.
% 
% Peaks are calculated after bandpass filtering df/f. Everything else is
% calculated after lowpass filtering df/f.
% 
% See source code for default values.
% Units for time and frequency are seconds and hertz respectively.
% 
% Normalization recipes:
%     df/f is calculated as (f - f0) / f1, where f0 and f1 are computed for
%       each time point using a function over a (moving) window.
%     Recipe for df/f:
%         configuration.f0Function = @movmean;
%         configuration.f1Function = @movmean;
%     Recipe for z-score:
%         configuration.f0Function = @movmean;
%         configuration.f1Function = @movstd;
% 	If dffEpochs covers the entire data set (e.g. [-Inf, Inf]), df/f is
%   calculated using a moving window. If dffEpochs covers a portion of
%   the data set (e.g. [10, 100]), df/f is calculated using a single window
%   in such period.
% 	
% See FPAexamples

% 2019-02-01. Leonardo Molina.
% 2020-02-19. Last modified.
function results = FPA(time, signal, reference, configuration)
    results.warnings = {};
    if nargin < 4
        configuration = struct();
    end
    
    % Read input configuration. Use defaults for missing parameters.
    configuration = setDefault(configuration, 'conditionEpochs', {'Condition A', [-Inf, Inf]});
    configuration = setDefault(configuration, 'dffEpochs', [-Inf, Inf]);
    configuration = setDefault(configuration, 'artifactEpochs', []);
    configuration = setDefault(configuration, 'bleachingEpochs', [-Inf, Inf]);
    configuration = setDefault(configuration, 'dffLowpassFrequency', 2.00);
    configuration = setDefault(configuration, 'peaksBandpassFrequency', [0.02, 0.20]);
    configuration = setDefault(configuration, 'bleachingLowpassFrequency', 0.1);
    configuration = setDefault(configuration, 'movingWindow', false);
    configuration = setDefault(configuration, 'f0Function', @movmean);
    configuration = setDefault(configuration, 'f0Window', 600);
    configuration = setDefault(configuration, 'f1Function', @movstd);
    configuration = setDefault(configuration, 'f1Window', 600);
    configuration = setDefault(configuration, 'thresholdingFunction', @mad);
    configuration = setDefault(configuration, 'thresholdFactor', 0.10);
    configuration = setDefault(configuration, 'triggeredWindow', 10);
    configuration = setDefault(configuration, 'f0', []);
    configuration = setDefault(configuration, 'f1', []);
    
    % Sampling frequency.
    sourceFrequency = 1 / mean(diff(time));
    configuration = setDefault(configuration, 'resamplingFrequency', sourceFrequency);
    
    % Make vectors of equal size.
    if isempty(reference)
        k = isnan(time) | isnan(signal);
        time(k) = [];
        signal(k) = [];
    else
        k = isnan(time) | isnan(signal) | isnan(reference);
        time(k) = [];
        signal(k) = [];
        reference(k) = [];
    end
    
    % Resample to target frequency.
    if configuration.resamplingFrequency < sourceFrequency
        % Express frequency as a ratio p/q.
        [p, q] = rat(configuration.resamplingFrequency / sourceFrequency);
        % Resample: interpolate every p/q/f, upsample by p, filter, downsample by q.
        [signal, time2] = resample(signal, time, configuration.resamplingFrequency, p, q);
        if ~isempty(reference)
            reference = resample(reference, time, configuration.resamplingFrequency, p, q);
        end
        time = time2;
    else
        results.warnings{end + 1} = warn('Cannot resample to frequencies higher than the source frequency (%.2f Hz).', sourceFrequency);
    end
    nSamples = numel(time);
    if isempty(reference)
        reference = ones(nSamples, 1);
    end
    
    % Replace artifacts with straight lines.
    % Index of all points.
    allIds = colon(1, nSamples)';
    % Index of artifacts and non-artifacts.
    badId = time2id(time, configuration.artifactEpochs);
    goodId = setdiff(allIds, badId);
    % Interpolate.
    reference2 = reference;
    signal2 = signal;
    reference2(badId) = interp1(goodId, reference(goodId), badId);
    signal2(badId) = interp1(goodId, signal(goodId), badId);
    
    % Remove high-frequency oscillations to detect bleaching decay.
    if configuration.bleachingLowpassFrequency <= configuration.resamplingFrequency / 2
        bleachingFilter = designfilt('lowpassiir', 'HalfPowerFrequency', configuration.bleachingLowpassFrequency, 'SampleRate', configuration.resamplingFrequency, 'DesignMethod', 'butter', 'FilterOrder', 12);
        rLowpass = filtfilt(bleachingFilter, reference2);
        sLowpass = filtfilt(bleachingFilter, signal2);
    else
        rLowpass = reference2;
        sLowpass = signal2;
        results.warnings{end + 1} = warn('[bleaching-correction] Cannot lowpass to frequencies higher than half of the resampling frequency (%.2f Hz).', configuration.resamplingFrequency / 2);
    end
    
    % Fit exponential decay at given epochs.
    bleachingCorrectionId = time2id(time, configuration.bleachingEpochs);
    rFit = fit(time(bleachingCorrectionId), rLowpass(bleachingCorrectionId), fittype('exp1'));
    sFit = fit(time(bleachingCorrectionId), sLowpass(bleachingCorrectionId), fittype('exp1'));
    rBleaching = rFit(time);
    sBleaching = sFit(time);
    
    % Correct for bleaching artifact.
    rCorrected = reference2 - rBleaching;
    sCorrected = signal2 - sBleaching;
    
    % Correct for movement artifacts.
    f = sCorrected - rCorrected;
    
    % Calculate df/f.
    if isempty(configuration.f0) && isempty(configuration.f1)
        dffIds = time2id(time, configuration.dffEpochs);
        if configuration.movingWindow
            % Compute individual baseline points for all data points using a moving window from all data points.
            f0 = configuration.f0Function(f, nanmin(round(configuration.f0Window * configuration.resamplingFrequency), numel(f)));
            f1 = configuration.f1Function(f, nanmin(round(configuration.f1Window * configuration.resamplingFrequency), numel(f)));
        else
            % Compute global baseline points from data the given epochs.
            f0 = configuration.f0Function(f(dffIds), numel(dffIds), 'Endpoints', 'discard');
            f1 = configuration.f1Function(f(dffIds), numel(dffIds), 'Endpoints', 'discard');
        end
    else
        f0 = configuration.f0;
        f1 = configuration.f1;
        if isempty(f0)
            f1 = f0;
        elseif isempty(f1)
            f0 = f1;
        end
    end
    dff = (f - f0) ./ f1;
    
    % Band-pass filter to detect peaks.
    if all(configuration.peaksBandpassFrequency <= configuration.resamplingFrequency / 2)
        bandpassFilter = designfilt('bandpassfir', 'CutoffFrequency1', configuration.peaksBandpassFrequency(1), 'CutoffFrequency2', configuration.peaksBandpassFrequency(2), 'SampleRate', configuration.resamplingFrequency, 'DesignMethod', 'window', 'FilterOrder', 12);
        dffBandpass = filtfilt(bandpassFilter, dff);
    else
        results.warnings{end + 1} = warn('[peak detection] Cannot bandpass at the given frequency range because at least one of the frequencies is greater than the resampling frequency (%.2f Hz).', configuration.resamplingFrequency);
        dffBandpass = dff;
    end
    
    % Get peak threshold.
    boundaryWindow =  ceil(configuration.peaksBandpassFrequency(1) * configuration.resamplingFrequency);
    useIds = setdiff(time2id(time, cat(2, configuration.conditionEpochs{2:2:end})), [1:boundaryWindow, numel(time) - boundaryWindow + 1]');
    peakThreshold = mean(dffBandpass(useIds)) + configuration.thresholdFactor * configuration.thresholdingFunction(dffBandpass(useIds));
    [~, peaksId] = findpeaks(+dffBandpass, 'MinPeakHeight', peakThreshold);
    peaksId = intersect(peaksId, useIds);
    valleyThreshold = -peakThreshold;
    [~, valleysId] = findpeaks(-dffBandpass, 'MinPeakHeight', peakThreshold);
    valleysId = intersect(valleysId, useIds);
    
    % Low-pass.
    if configuration.dffLowpassFrequency <= configuration.resamplingFrequency / 2
        lowpassFilter = designfilt('lowpassiir', 'HalfPowerFrequency', configuration.dffLowpassFrequency, 'SampleRate', configuration.resamplingFrequency, 'DesignMethod', 'butter', 'FilterOrder', 12);
        dffLowpass = filtfilt(lowpassFilter, dff);
    else
        results.warnings{end + 1} = warn('[dff] Cannot lowpass to frequencies smaller than half of the resampling frequency (%.2f Hz).', configuration.resamplingFrequency / 2);
        dffLowpass = dff;
    end
    
    % Split peaks/traces by conditions.
    % Number of samples in a triggered window.
    window = round(configuration.triggeredWindow * configuration.resamplingFrequency);
    % Force odd count.
    if mod(window, 2) == 0
        window = window - 1;
    end
    % Index template to apply around each peak.
    windowTemplate = -(window - 1) / 2:(window - 1) / 2;
    % Filter out out-of-range traces.
    peaksId = peaksId(peaksId > (window - 1) / 2 & peaksId + (window - 1) / 2 < nSamples);
    nPeaks = numel(peaksId);
    % Index of all triggered traces; one per row.
    if nPeaks > 0
        triggeredId = bsxfun(@plus, windowTemplate, peaksId);
    else
        triggeredId = [];
    end
    % Assign a group to each trace (rows in triggeredId) according to the epoch they are.
    peakGroups = zeros(1, nPeaks);
    % Boolean index corresponding to each condition.
    nEpochs = numel(configuration.conditionEpochs) / 2;
    epochBool = cell(1, nEpochs);
    for e = 1:nEpochs
        % Epoch index.
        ids = time2id(time, configuration.conditionEpochs{2 * e});
        epochBool{e} = ismember(allIds, ids);
        % Assign group.
        k = epochBool{e}(peaksId);
        peakGroups(k) = e;
    end
    uniqueGroups = unique(peakGroups);
    % An epoch may not have any peaks.
    uniqueGroups = uniqueGroups(uniqueGroups > 0);
    nGroups = numel(uniqueGroups);
    
    % Style.
    cmap = lines();
    xlims = time([1, end]);
    
    % Plot raw signal and bleaching.
    results.figures = figure('name', 'FPA: df/f');
    ax.raw = subplot(3, 1, 1);
    ax.raw.XTick = [];
    hold(ax.raw, 'all');
    yy = [signal(:); reference(:); sBleaching(:)];
    ylims = [min(yy), max(yy)];
    plotEpochs(configuration.conditionEpochs, xlims, ylims, cmap, true);
    plot(ax.raw, time, signal, 'DisplayName', 'Signal');
    plot(ax.raw, time, reference, 'DisplayName', 'Reference');
    plot(ax.raw, time, sBleaching, 'Color', [0, 0, 0], 'LineStyle', '--', 'DisplayName', 'Bleaching');
    axis(ax.raw, 'tight');
    legend(ax.raw, 'show');
    
    % Plot band-pass filtered signal and peaks.
    ax.peaks = subplot(3, 1, 2);
    ax.peaks.XTick = [];
    hold(ax.peaks, 'all');
    yy = dffBandpass(:);
    ylims = [min(yy), max(yy)]; %[prctile(yy, 1), max(peakThreshold, prctile(yy, 99))];
    ylims(1) = ylims(1) - 0.25 * diff(ylims);
    ylims(2) = ylims(2) + 0.25 * diff(ylims);
    plotEpochs(configuration.conditionEpochs, xlims, ylims, cmap, false);
    plot(ax.peaks, time, dffBandpass, 'Color', [0, 0, 0], 'DisplayName', 'band-pass df/f');
    plot(ax.peaks, time(peaksId), dffBandpass(peaksId), 'Color', [1, 0, 0], 'LineStyle', 'none', 'Marker', 'o', 'DisplayName', sprintf('%i peaks', numel(peaksId)));
    plot(ax.peaks, time([1, end]), peakThreshold([1, 1]), 'Color', [0, 0, 0], 'LineStyle', '--', 'DisplayName', 'threshold');
    plot(ax.peaks, time(valleysId), dffBandpass(valleysId), 'Color', [0, 0, 0], 'LineStyle', 'none', 'Marker', 'o', 'DisplayName', sprintf('%i valleys', numel(valleysId)));
    plot(ax.peaks, time([1, end]), valleyThreshold([1, 1]), 'Color', [0, 0, 0], 'LineStyle', '--', 'HandleVisibility', 'off');
    ylim(ylims);
    legend(ax.peaks, 'show');
    
    % Plot corrected signal, low-pass filtered signal and peaks.
    ax.processed = subplot(3, 1, 3);
    hold(ax.processed, 'all');
    yy = [dff(:); dffLowpass(:)];
    ylims = [min(yy), max(yy)]; %[prctile(yy, 1), max(peakThreshold, prctile(yy, 99))];
    ylims(1) = ylims(1) - 0.25 * diff(ylims);
    ylims(2) = ylims(2) + 0.25 * diff(ylims);
    plotEpochs(configuration.conditionEpochs, xlims, ylims, cmap, false);
    plot(ax.processed, time, dff, 'DisplayName', 'df/f');
    plot(ax.processed, time, dffLowpass, 'Color', [0, 0, 0], 'DisplayName', 'low-pass df/f');
    plot(ax.processed, time(peaksId), dffLowpass(peaksId), 'Color', [1, 0, 0], 'LineStyle', 'none', 'Marker', 'o', 'DisplayName', 'peaks');
    plot(ax.processed, time(valleysId), dffLowpass(valleysId), 'Color', [0, 0, 0], 'LineStyle', 'none', 'Marker', 'o', 'DisplayName', 'valleys');
    ylim(ylims);
    legend(ax.processed, 'show');
    
    % Move axes together.
    linkaxes([ax.raw, ax.peaks, ax.processed], 'x');
    xlim(ax.raw, [time(1), time(end)]);
    
    xlabel('Time (s)');
    ylabel('df/f');
    
    % Plot power spectrum.
    results.figures(end + 1) = figure('name', 'FPA: Power spectrum');
    axs = cell(1, nEpochs);
    for e = 1:nEpochs
        axs{e} = subplot(nEpochs, 1, e);
        epochName = configuration.conditionEpochs{2 * e - 1};
        ids = time2id(time, configuration.conditionEpochs{2 * e});
        d = dff(ids);
        n = numel(ids);
        halfN = floor(n / 2);
        f = fft(d);
        % Two-sided spectrum.
        p2 = abs(f / n);
        % Single-sided amplitude spectrum.
        p1 = p2(1:halfN + 1);
        p1(2:end-1) = 2 * p1(2:end-1);
        % Create frequency vector for range.
        fs = configuration.resamplingFrequency * (0:halfN) / n;
        plot(fs, p1);
        title(sprintf('%s - Power spectrum', epochName));
    end
    ylabel('Power');
    xlabel('Frequency (Hz)');
    linkaxes([axs{:}], 'x');
    
    % Plot triggered average.
    results.figures(end + 1) = figure('name', 'FPA: Triggered average');
    ax.trigger = axes();
    hold(ax.trigger, 'all');
    timeTemplate = windowTemplate / configuration.resamplingFrequency;
    for e = 1:nGroups
        group = uniqueGroups(e);
        triggeredDff = dff(triggeredId(peakGroups == group, :));
        triggeredDff = reshape(triggeredDff, numel(triggeredDff) / window, window);
        triggeredMean = mean(triggeredDff, 1);
        h1 = plot(timeTemplate, triggeredMean, 'HandleVisibility', 'off');
        epochName = configuration.conditionEpochs{2 * e - 1};
        triggeredSem = std(triggeredDff, [], 1) / sqrt(size(triggeredDff, 1));
        semAtZero = triggeredSem(ceil(window / 2));
        nPeaks = sum(peakGroups == group);
        text = sprintf('%s (SEM=%.4f, n = %i)', epochName, semAtZero, nPeaks);
        vertices = [timeTemplate; triggeredMean + triggeredSem / 2];
        vertices = cat(2, vertices, [fliplr(timeTemplate); fliplr(triggeredMean - triggeredSem / 2)])';
        faces = 1:2 * window;
        patch('Faces', faces, 'Vertices', vertices, 'FaceColor', h1.Color, 'EdgeColor', 'none', 'FaceAlpha', 0.10, 'DisplayName', text);
    end
    title('Triggered average');
    legend('show');
    xlabel('Time (s)');
    ylabel('df/f');
    axis(ax.trigger, 'tight');
    
    % Boxplot of dff.
    results.figures(end + 1) = figure('name', 'FPA: Boxplot');
    epochIds = zeros(0, 1);
    epochGroups = zeros(0, 1);
    epochLabels = cell(1, nEpochs);
    for e = 1:nEpochs
        ids = time2id(time, configuration.conditionEpochs{2 * e});
        epochIds = cat(1, epochIds, ids);
        thisEpochGroups = repmat(e, [numel(ids), 1]);
        epochGroups = cat(1, epochGroups, thisEpochGroups);
        epochLabels{e} = sprintf('%s (STD:%.4f)', configuration.conditionEpochs{2 * e - 1}, std(dff(ids)));
    end
    boxplot(dff(epochIds), epochGroups, 'Labels', epochLabels);
    title('Stats on df/f traces for each condition');
    ylabel('df/f');
    
    results.time = time;
    results.dff = dff;
    results.peaksId = peaksId;
    results.reference = reference2;
    results.signal = signal2;
    results.epochIds = epochIds;
    results.epochGroups = epochGroups;
    results.f0 = mean(f0);
    results.f1 = mean(f1);
end

function configuration = setDefault(configuration, fieldname, value)
    if ~isfield(configuration, fieldname)
        configuration.(fieldname) = value;
    end
end

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

function output = warn(format, varargin)
    output = sprintf('[%s] %s', mfilename(), sprintf(format, varargin{:}));
end