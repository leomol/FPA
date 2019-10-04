% FPA(time, signal, reference, configuration)
% 
% Correct signal from bleaching and artifacts, normalize and detect peaks of
% spontaneous activity based on the given parameters. Signal and reference
% are column vectors.
% 
% Overall analysis steps:
% 	-Resample data.
% 	-Fit an exponential decay in a portion of the data representative of bleaching.
% 	-Normalize signal to reference.
% 	-Compute df/f in a (moving) window.
% 	-Find peaks of spontaneous activity.
% 	-Plot corrected signal and peaks; highlight epochs.
% 	-Plot triggered averages.
% 	-Plot power spectrum.
% 	-Plot stats.
% 
% configuration is a struct with the following fields (defaults are used for missing fields):
%     conditionEpochs - Epochs for different conditions: {'epoch1', [start1, end1, start2, end2, ...], 'epoch2', ...}
%     baselineEpochs - Time epochs (s) to include for df/f normalization.
%     bleachingEpochs - Time epochs (s) to include for bleaching correction.
%     artifactEpochs - Time epochs (s) to remove.
%     resamplingFrequency - Resampling frequency (Hz).
%     dffLowpassFrequency - Lowpass frequency to filter df/f.
%     peaksBandpassFrequency - Low/High frequencies to compute peaks.
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
% 	If baselineEpochs covers the entire data set (e.g. [-Inf, Inf]), df/f is
%   calculated using a moving window. If baselineEpochs covers a portion of
%   the data set (e.g. [10, 100]), df/f is calculated using a single window
%   in such period.
% 	
% See FPAexamples

% 2019-02-01. Leonardo Molina.
% 2019-10-04. Last modified.
function results = FPA(time, signal, reference, configuration)
    if nargin < 4
        configuration = struct();
    end
    
    % Read input configuration. Use defaults for missing parameters.
    configuration = setDefault(configuration, 'conditionEpochs', {'Condition A', [-Inf, Inf]});
    configuration = setDefault(configuration, 'baselineEpochs', [-Inf, Inf]);
    configuration = setDefault(configuration, 'artifactEpochs', []);
    configuration = setDefault(configuration, 'bleachingEpochs', [-Inf, Inf]);
    configuration = setDefault(configuration, 'resamplingFrequency', 50);
    configuration = setDefault(configuration, 'dffLowpassFrequency', 0.2);
    configuration = setDefault(configuration, 'peaksBandpassFrequency', [0.02, 0.2]);
    configuration = setDefault(configuration, 'bleachingLowpassFrequency', 0.1);
    configuration = setDefault(configuration, 'f0Function', @movmean);
    configuration = setDefault(configuration, 'f0Window', 10);
    configuration = setDefault(configuration, 'f1Function', @movstd);
    configuration = setDefault(configuration, 'f1Window', 10);
    configuration = setDefault(configuration, 'thresholdingFunction', @mad);
    configuration = setDefault(configuration, 'thresholdFactor', 0.10);
    configuration = setDefault(configuration, 'triggeredWindow', 10);
    
    % Sampling frequency.
    sourceFrequency = 1 / mean(diff(time));
    
    % Resample to target frequency.
    % Express frequency as a ratio p/q.
    [p, q] = rat(configuration.resamplingFrequency / sourceFrequency);
    % Resample: interpolate every p/q/f, upsample by p, filter, downsample by q.
    reference = resample(reference, time, configuration.resamplingFrequency, p, q);
    [signal, time] = resample(signal, time, configuration.resamplingFrequency, p, q);
    nSamples = numel(time);
    
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
    bleachingFilter = designfilt('lowpassiir', 'HalfPowerFrequency', configuration.bleachingLowpassFrequency, 'SampleRate', configuration.resamplingFrequency, 'DesignMethod', 'butter', 'FilterOrder', 12);
    rLowpass = filtfilt(bleachingFilter, reference2);
    sLowpass = filtfilt(bleachingFilter, signal2);
    
    % Fit exponential decay at given epochs.
    bleachingCorrectionId = time2id(time, configuration.bleachingEpochs);
    rFit = fit(time(bleachingCorrectionId), rLowpass(bleachingCorrectionId), fittype('exp1'));
    sFit = fit(time(bleachingCorrectionId), sLowpass(bleachingCorrectionId), fittype('exp1'));
    rBleaching = rFit(time);
    sBleaching = sFit(time);
    
    % Correct for bleaching artifact.
    rCorrected = reference2 ./ rBleaching;
    sCorrected = signal2 ./ sBleaching;
    
    % Correct for movement artifacts.
    f = sCorrected ./ rCorrected;
    
    % Calculate df/f.
    baselineId = time2id(time, configuration.baselineEpochs);
    if all(ismember(allIds, baselineId))
        % Compute baseline from a moving window.
        f0 = configuration.f0Function(f, nanmin(round(configuration.f0Window * configuration.resamplingFrequency), numel(f)));
        f1 = configuration.f1Function(f, nanmin(round(configuration.f1Window * configuration.resamplingFrequency), numel(f)));
    else
        % Compute baseline at given epoch.
        f0 = configuration.f0Function(f(baselineId), numel(baselineId), 'Endpoints', 'discard');
        f1 = configuration.f1Function(f(baselineId), numel(baselineId), 'Endpoints', 'discard');
    end
    dff = (f - f0) ./ f1;
    
    % Band-pass filter to detect peaks.
    bandpassFilter = designfilt('bandpassiir', 'HalfPowerFrequency1', configuration.peaksBandpassFrequency(1), 'HalfPowerFrequency2', configuration.peaksBandpassFrequency(2), 'SampleRate', configuration.resamplingFrequency, 'DesignMethod', 'butter', 'FilterOrder', 12);
    dffBandpass = filtfilt(bandpassFilter, dff);
    
    % Get peak threshold.
    boundaryWindow =  ceil(configuration.peaksBandpassFrequency(1) * configuration.resamplingFrequency);
    cleanIds = allIds(boundaryWindow:end - boundaryWindow + 1);
    peakThreshold = mean(dffBandpass(cleanIds)) + configuration.thresholdFactor * configuration.thresholdingFunction(dffBandpass(cleanIds));
    [~, peaksId] = findpeaks(dffBandpass, 'MinPeakHeight', peakThreshold);
    
    % Low-pass.
    lowpassFilter = designfilt('lowpassiir', 'HalfPowerFrequency', configuration.dffLowpassFrequency, 'SampleRate', configuration.resamplingFrequency, 'DesignMethod', 'butter', 'FilterOrder', 12);
    dffLowpass = filtfilt(lowpassFilter, dff);
    
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
        epochIds = time2id(time, configuration.conditionEpochs{2 * e});
        epochBool{e} = ismember(allIds, epochIds);
        % Assign group.
        k = epochBool{e}(peaksId);
        peakGroups(k) = e;
    end
    uniqueGroups = unique(peakGroups);
    % An epoch may not have any peaks.
    uniqueGroups = uniqueGroups(uniqueGroups > 0);
    nGroups = numel(uniqueGroups);
    
    % Color palette.
    cmap = lines();
    
    % Plot raw signal and bleaching.
    figure('name', 'FPA: df/f');
    ax.raw = subplot(3, 1, 1);
    ax.raw.XTick = [];
    hold(ax.raw, 'all');
    yy = [signal(:); reference(:); sBleaching(:)];
    ylims = [min(yy), max(yy)];
    plotEpochs(configuration.conditionEpochs, ylims, cmap, true);
    plot(ax.raw, time, signal, 'DisplayName', 'Signal');
    plot(ax.raw, time, reference, 'DisplayName', 'Reference');
    plot(ax.raw, time, sBleaching, 'Color', [0, 0, 0], 'LineStyle', '--', 'DisplayName', 'Bleaching');
    legend(ax.raw, 'show');
    
    % Plot band-pass filtered signal and peaks.
    ax.peaks = subplot(3, 1, 2);
    ax.peaks.XTick = [];
    hold(ax.peaks, 'all');
    yy = dffBandpass(:);
    ylims = [prctile(yy, 5), max(peakThreshold, prctile(yy, 95))];
    ylims(1) = ylims(1) - 0.25 * diff(ylims);
    ylims(2) = ylims(2) + 0.25 * diff(ylims);
    plotEpochs(configuration.conditionEpochs, ylims, cmap, false);
    plot(ax.peaks, time, dffBandpass, 'Color', [0, 0, 0], 'DisplayName', 'band-pass df/f');
    plot(ax.peaks, time(peaksId), dffBandpass(peaksId), 'ro', 'DisplayName', 'peaks');
    plot(ax.peaks, time([1, end]), peakThreshold([1, 1]), 'Color', [0, 0, 0], 'LineStyle', '--', 'DisplayName', 'threshold');
    ylim(ylims);
    legend(ax.peaks, 'show');
    
    % Plot corrected signal, low-pass filtered signal and peaks.
    ax.processed = subplot(3, 1, 3);
    ax.processed.XTick = [];
    hold(ax.processed, 'all');
    yy = [dff(:); dffLowpass(:)];
    ylims = [min(yy), max(yy)];
    plotEpochs(configuration.conditionEpochs, ylims, cmap, false);
    plot(ax.processed, time, dff, 'DisplayName', 'df/f');
    plot(ax.processed, time, dffLowpass, 'Color', [0, 0, 0], 'DisplayName', 'low-pass df/f');
    plot(ax.processed, time(peaksId), dffLowpass(peaksId), 'Color', [1, 0, 0], 'LineStyle', 'none', 'Marker', 'o', 'DisplayName', 'peaks');
    legend(ax.processed, 'show');
    
    % Move axes together and use all space.
    linkaxes([ax.raw, ax.peaks, ax.processed], 'x');
    axis([ax.raw, ax.processed], 'tight');
    
    xlabel('Time (s)');
    ylabel('df/f');
    
    % Plot power spectrum.
    figure('name', 'FPA: Power spectrum');
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
    figure('name', 'FPA: Triggered average');
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
    figure('name', 'FPA: Boxplot');
    allEpochIds = zeros(0, 1);
    allEpochGroups = zeros(0, 1);
    epochLabels = cell(1, nEpochs);
    for e = 1:nEpochs
        epochIds = time2id(time, configuration.conditionEpochs{2 * e});
        allEpochIds = cat(1, allEpochIds, epochIds);
        epochGroups = repmat(e, [numel(epochIds), 1]);
        allEpochGroups = cat(1, allEpochGroups, epochGroups);
        epochLabels{e} = sprintf('%s (STD:%.4f)', configuration.conditionEpochs{2 * e - 1}, std(dff(epochIds)));
    end
    boxplot(dff(allEpochIds), allEpochGroups, 'Labels', epochLabels);
    title('Stats on df/f traces for each condition');
    ylabel('df/f');
    
    results.time = time;
    results.dff = dff;
    results.peaksId = peaksId;
    results.reference = reference2;
    results.signal = signal2;
end

function configuration = setDefault(configuration, fieldname, value)
    if ~isfield(configuration, fieldname)
        configuration.(fieldname) = value;
    end
end

function plotEpochs(epochs, ylims, cmap, show)
    for e = 1:numel(epochs) / 2
        epochName = epochs{2 * e - 1};
        [faces, vertices] = patchEpochs(epochs{2 * e}, ylims(1), ylims(2));
        if show
            patch('Faces', faces, 'Vertices', vertices, 'FaceColor', cmap(e, :), 'EdgeColor', 'none', 'FaceAlpha', 0.50, 'DisplayName', sprintf('%s', epochName));
        else
            patch('Faces', faces, 'Vertices', vertices, 'FaceColor', cmap(e, :), 'EdgeColor', 'none', 'FaceAlpha', 0.50, 'HandleVisibility', 'off');
        end
    end
end