% FPA(time, signal, reference, configuration)
% 
% Correct signal from bleaching and artifacts; detect peaks of spontaneous
% activity based on the given parameters.
% Signal and reference are column vectors.
% 
% Overall analysis steps:
%   -Resample to target frequency.
%   -Remove artifacts: Replace artifacts with lines in flagged regions.
% 	-Correct for bleaching: Fit an exponential decay in low-pass data.
% 	-Correct for motion artifacts: Subtract bleaching corrected signals.
% 	-Compute df/f or z-score according to settings.
% 	-Find peaks of spontaneous activity in low-pass signal.
%   -Plot 1:
%     -Raw signal and bleaching fit.
%     -Low-pass signal and peaks.
%     -Motion and bleaching corrected signal.
%   -Plot 2:
%     -Power spectrum of each epoch.
%   -Plot 3:
%     -Boxplot.
%   -Plot 4:
%     -Triggered average of spontaneous activity (if any peaks are found).
% 
% configuration is a struct with the following fields (defaults are used for missing fields):
%     conditionEpochs - Epochs for different conditions: {'epoch1', [start1, end1, start2, end2, ...], 'epoch2', ...}
%     bleachingEpochs - Time epochs (s) to include for bleaching correction.
%     artifactEpochs - Time epochs (s) to remove.
%     resamplingFrequency - Resampling frequency (Hz).
%     dffLowpassFrequency - Lowpass frequency to filter df/f.
%     peaksLowpassFrequency - Low/High frequencies to detect peaks.
%     bleachingLowpassFrequency - Lowpass frequency to detect bleaching decay.
%     thresholdingFunction - One of @mad, @std.
%     thresholdFactor - Thresholding cut-off.
%     triggeredWindow - Length of time to capture around each peak of spontaneous activity.
%     dffEpochs - Time epochs (s) to include for df/f normalization.
%     f0Function - One of @movmean, @movmedian, @movmin.
%     f0Window - Length of the moving window to calculate f0.
%     f1Function - One of @movmean, @movmedian, @movmin, @movstd.
%     f1Window - Length of the moving window to calculate f1.
% 
% See source code for detailed analysis steps and default parameters.
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
%     
%     Example 1:
%         Use the first minute as the baseline for the rest of the data.
%         configuration.f0Function = @movmean;
%         configuration.f1Function = @movstd;
%         configuration.dffEpochs = [0, 60]
%         >>> f0Window and f1Window are ignored because dffEpochs takes precedence.
%     Example 2:
%         Compute a moving baseline that is 1min wide.
%         configuration.f0Function = @movmean;
%         configuration.f1Function = @movstd;
%         configuration.f0Window = 60;
%         configuration.f1Window = 60;
%         >>> dffEpochs must be empty (configuration.dffEpochs = [];) which is the default.
%     Example 3:
%         Compute a common baseline from the entire recording.
%         configuration.f0Function = @movmean;
%         configuration.f1Function = @movstd;
%         configuration.f0Window = Inf;
%         configuration.f1Window = Inf;
%         >>> dffEpochs must be empty.
% 	
% See examples

% 2019-02-01. Leonardo Molina.
% 2020-10-29. Last modified.
function results = FPA(time, signal, reference, configuration)
    results.warnings = {};
    if nargin < 4
        configuration = struct();
    end
    
    % Read input configuration. Use defaults for missing parameters.
    configuration = setDefault(configuration, 'conditionEpochs', {'Condition A', [-Inf, Inf]});
    configuration = setDefault(configuration, 'dffEpochs', []);
    configuration = setDefault(configuration, 'bleachingEpochs', [-Inf, Inf]);
    configuration = setDefault(configuration, 'artifactEpochs', []);
    configuration = setDefault(configuration, 'dffLowpassFrequency', 2.0);
    configuration = setDefault(configuration, 'peaksLowpassFrequency', 0.2);
    configuration = setDefault(configuration, 'bleachingLowpassFrequency', 0.1);
    configuration = setDefault(configuration, 'f0Function', @movmean);
    configuration = setDefault(configuration, 'f0Window', 600);
    configuration = setDefault(configuration, 'f1Function', @movstd);
    configuration = setDefault(configuration, 'f1Window', 600);
    configuration = setDefault(configuration, 'thresholdingFunction', @mad);
    configuration = setDefault(configuration, 'thresholdFactor', 0.10);
    configuration = setDefault(configuration, 'triggeredWindow', 10);
    configuration = setDefault(configuration, 'f0', []);
    configuration = setDefault(configuration, 'f1', []);
    configuration = setDefault(configuration, 'plot', true);
    
    percentile = 0.99;
    grow = 0.50;
    
    % Sampling frequency.
    sourceFrequency = 1 / median(diff(time));
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
    elseif configuration.resamplingFrequency ~= sourceFrequency
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
        if isempty(dffIds)
            % For each data point, compute a baseline from its neighbors.
            n0 = nanmin(round(configuration.f0Window * configuration.resamplingFrequency), numel(f));
            f0 = configuration.f0Function(f, n0);
            n1 = nanmin(round(configuration.f1Window * configuration.resamplingFrequency), numel(f));
            f1 = configuration.f1Function(f, n1);
        else
            % For all points, compute a common baseline from the given epochs.
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
    
    % Low-pass filter to detect peaks.
    if configuration.peaksLowpassFrequency <= configuration.resamplingFrequency / 2
        peaksFilter = designfilt('lowpassiir', 'HalfPowerFrequency', configuration.peaksLowpassFrequency, 'SampleRate', configuration.resamplingFrequency, 'DesignMethod', 'butter', 'FilterOrder', 12);
        peaksLowpass = filtfilt(peaksFilter, dff);
    else
        results.warnings{end + 1} = warn('[peak detection] Cannot lowpass to frequencies smaller than half of the resampling frequency (%.2f Hz).', configuration.resamplingFrequency / 2);
        peaksLowpass = dff;
    end
    
    % Get peak threshold.
    excludeWindow =  ceil(configuration.peaksLowpassFrequency * configuration.resamplingFrequency);
    useIds = setdiff(time2id(time, cat(2, configuration.conditionEpochs{2:2:end})), [1:excludeWindow, numel(time) - excludeWindow + 1]');
    peakThreshold = mean(peaksLowpass(useIds)) + configuration.thresholdFactor * configuration.thresholdingFunction(peaksLowpass(useIds));
    valleyThreshold = -peakThreshold;
    if any(peaksLowpass >= peakThreshold)
        [~, peaksId] = findpeaks(+peaksLowpass, 'MinPeakHeight', peakThreshold);
        peaksId = intersect(peaksId, useIds);
    else
        peaksId = [];
    end
    if any(-peaksLowpass >= peakThreshold)
        [~, valleysId] = findpeaks(-peaksLowpass, 'MinPeakHeight', peakThreshold);
        valleysId = intersect(valleysId, useIds);
    else
        valleysId = [];
    end
    
    % Low-pass.
    if configuration.dffLowpassFrequency <= configuration.resamplingFrequency / 2
        lowpassFilter = designfilt('lowpassiir', 'HalfPowerFrequency', configuration.dffLowpassFrequency, 'SampleRate', configuration.resamplingFrequency, 'DesignMethod', 'butter', 'FilterOrder', 12);
        dffLowpass = filtfilt(lowpassFilter, dff);
    else
        results.warnings{end + 1} = warn('[dff lowpass] Cannot lowpass to frequencies smaller than half of the resampling frequency (%.2f Hz).', configuration.resamplingFrequency / 2);
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
    % Epochs may not have peaks.
    uniqueGroups = uniqueGroups(uniqueGroups > 0);
    nGroups = numel(uniqueGroups);
    
    epochIds = zeros(0, 1);
    epochGroups = zeros(0, 1);
    epochStatLabels = cell(1, nEpochs);
    area = zeros(1, nEpochs);
    peaksCount = zeros(1, nEpochs);
    valleysCount = zeros(1, nEpochs);
    for e = 1:nEpochs
        ids = time2id(time, configuration.conditionEpochs{2 * e});
        area(e) = sum(dff(ids));
        peaksCount(e) = sum(ismember(ids, peaksId));
        valleysCount(e) = sum(ismember(ids, valleysId));
        epochIds = cat(1, epochIds, ids);
        thisEpochGroups = repmat(e, [numel(ids), 1]);
        epochGroups = cat(1, epochGroups, thisEpochGroups);
        epochStatLabels{e} = sprintf('\nmean:%.2f\nstd:%.2f', mean(dff(ids)), std(dff(ids)));
    end
    
    % Style.
    cmap = lines();
    xlims = time([1, end]);
    
    if configuration.plot
        signalColor = [0 0.4470 0.7410];
        referenceColor = [0.8500 0.3250 0.0980];
        lowpassColor = [0.4940 0.1840 0.5560];
        peaksLineColor = [0.4660 0.6740 0.1880];
        peaksMarkerColor = [1, 0, 0];
        dashColor = [0, 0, 0];
        
        results.figures = [];
        
        results.figures(end + 1) = figure('name', 'FPA: df/f');
        
        % Plot signal/reference resampled and bleaching model.
        ax.raw = subplot(4, 1, 1);
        ax.raw.XTick = [];
        hold(ax.raw, 'all');
        yy = [signal(:); reference(:); sBleaching(:)];
        ylims = limits(yy, percentile, grow);
        plotEpochs(configuration.conditionEpochs, xlims, ylims, cmap, true);
        plot(ax.raw, time, signal, 'Color', signalColor, 'DisplayName', 'Signal');
        plot(ax.raw, time, reference, 'Color', referenceColor, 'DisplayName', 'Reference');
        plot(ax.raw, time, sBleaching, 'Color', dashColor, 'LineStyle', '--', 'DisplayName', 'Bleaching');
        ylim(ylims);
        legend(ax.raw, 'show');
        
        % Plot signal/reference resampled and corrected for bleaching.
        ax.corrected = subplot(4, 1, 2);
        ax.corrected.XTick = [];
        hold(ax.corrected, 'all');
        yy = [sCorrected(:); rCorrected(:);];
        ylims = limits(yy, percentile, grow);
        plotEpochs(configuration.conditionEpochs, xlims, ylims, cmap, false);
        plot(ax.corrected, time, sCorrected, 'Color', signalColor, 'DisplayName', 'Signal');
        plot(ax.corrected, time, rCorrected, 'Color', referenceColor, 'DisplayName', 'Reference');
        ylim(ylims);
        legend(ax.corrected, 'show');
        
        % Plot df/f.
        ax.filtered = subplot(4, 1, 3);
        ax.filtered.XTick = [];
        hold(ax.filtered, 'all');
        yy = [dff(:); peaksLowpass(:); dffLowpass(:)];
        ylims = limits(yy, percentile, grow);
        epochs = configuration.conditionEpochs;
        epochs(1:2:end) = arrayfun(@(e) sprintf('area:%.2f', area(e)), 1:nEpochs, 'UniformOutput', false);
        plotEpochs(epochs, xlims, ylims, cmap, true);
        plot(ax.filtered, time, dff, 'Color', signalColor, 'DisplayName', 'df/f');
        plot(ax.filtered, time, dffLowpass, 'Color', lowpassColor, 'DisplayName', sprintf('Lowpass filtered df/f (%.2fHz)', configuration.dffLowpassFrequency));
        plot(ax.filtered, time, peaksLowpass, 'Color', peaksLineColor, 'DisplayName', sprintf('Lowpass filtered df/f (%.2fHz)', configuration.peaksLowpassFrequency));
        ylim(ylims);
        legend(ax.filtered, 'show');

        % Plot df/f and peaks.
        ax.processed = subplot(4, 1, 4);
        hold(ax.processed, 'all');
        yy = peaksLowpass(:);
        ylims = limits(yy, percentile, grow);
        
        epochs = configuration.conditionEpochs;
        epochs(1:2:end) = arrayfun(@(e) sprintf('%i peaks / %i valleys', peaksCount(e), valleysCount(e)), 1:nEpochs, 'UniformOutput', false);
        plotEpochs(epochs, xlims, ylims, cmap, true);
        plot(ax.processed, time, peaksLowpass, 'Color', peaksLineColor, 'DisplayName', sprintf('Lowpass filtered df/f (%.2fHz)', configuration.peaksLowpassFrequency));
        plot(ax.processed, time([1, end]), peakThreshold([1, 1]), 'Color', dashColor, 'LineStyle', '--', 'DisplayName', 'threshold');
        plot(ax.processed, time([1, end]), valleyThreshold([1, 1]), 'Color', dashColor, 'LineStyle', '--', 'HandleVisibility', 'off');
        plot(ax.processed, time(peaksId), peaksLowpass(peaksId), 'Color', peaksMarkerColor, 'LineStyle', 'none', 'Marker', 'o', 'HandleVisibility', 'off');
        plot(ax.processed, time(valleysId), peaksLowpass(valleysId), 'Color', peaksMarkerColor, 'LineStyle', 'none', 'Marker', 'o', 'HandleVisibility', 'off');
        ylim(ylims);
        legend(ax.processed, 'show');

        % Move axes together.
        linkaxes(findobj(gcf(), 'type', 'axes'), 'x');
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
            ylim(limits(p1, percentile, grow));
            title(sprintf('%s - Power spectrum', epochName));
        end
        ylabel('Power');
        xlabel('Frequency (Hz)');
        linkaxes(findobj(gcf(), 'type', 'axes'), 'x');

        % Boxplot of dff.
        results.figures(end + 1) = figure('name', 'FPA: Boxplot');
        h = axes();
        boxplot(dff(epochIds), epochGroups, 'Labels',  configuration.conditionEpochs(1:2:end));
        hold('all');
        area = zeros(1, nEpochs);
        ylims = ylim();
        for e = 1:nEpochs
            ids = time2id(time, configuration.conditionEpochs{2 * e});
            area(e) = mean(dff(ids));
            text(e, ylims(2), epochStatLabels{e}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'Top');
        end
        plot(h.XTick, area, 'k.', 'DisplayName', 'Mean');
        ylabel('df/f');
        xtickangle(45);
        title('Stats on df/f traces for each condition');

        % Plot triggered average.
        results.figures(end + 1) = figure('name', 'FPA: Triggered average');
        ax.trigger = axes();
        hold(ax.trigger, 'all');
        if nGroups > 0
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
                label = sprintf('%s (SEM=%.4f, n = %i)', epochName, semAtZero, nPeaks);
                vertices = [timeTemplate; triggeredMean + triggeredSem / 2];
                vertices = cat(2, vertices, [fliplr(timeTemplate); fliplr(triggeredMean - triggeredSem / 2)])';
                faces = 1:2 * window;
                patch('Faces', faces, 'Vertices', vertices, 'FaceColor', h1.Color, 'EdgeColor', 'none', 'FaceAlpha', 0.10, 'DisplayName', label);
            end
        else
            text(ax.trigger, 0.5, 0.5, sprintf('No peaks above threshold %.2f (factor:%.2f)', peakThreshold, configuration.thresholdFactor), 'HorizontalAlignment', 'center');
        end
        title('Triggered average');
        legend('show');
        xlabel('Time (s)');
        ylabel('df/f');
        axis(ax.trigger, 'tight');
    end

    results.time = time;
    results.dff = dff;
    results.dffLowpass = dffLowpass;
    results.peaksLowpass = peaksLowpass;
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

function ylims = limits(x, percentile, grow)
    ylims = [prctile(x, 100 * (1 - percentile)), prctile(x, 100 * percentile)];
    delta = diff(ylims) * grow;
    ylims = [ylims(1) - delta, ylims(2) + delta];
end