% results = FPA(time, signal, reference, configuration)
% 
% Correct data from photobleaching and motion artifacts; normalize, filter,
% and detect peaks of spontaneous activity in different epochs.
% 
% Time, signal, and reference are column vectors.
% 
% Processing steps:
%   -Resample data to target frequency.
%   -Replace artifacts with linear interpolation in flagged regions.
% 	-Correct for photobleaching by modeling an exponential decay in low-pass filtered data.
% 	-Correct for motion artifacts by subtracting reference to signal, after a polynomial fit.
%   -Remove fast oscillations with a low-pass filter.
% 	-Normalize data as df/f or z-score according to settings.
% 	-Find peaks of spontaneous activity in low-pass filtered data.
%   -Plot 1:
%     -Raw signal and reference, and photobleaching model.
%     -Signal and reference corrected for photobleaching.
%     -Baseline correction.
%     -Normalization.
%     -Peaks.
%   -Plot 2:
%     -Power spectrum for each epoch.
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
%     lowpassFrequency - Lowest frequency permitted in signal.
%     peaksLowpassFrequency - Lowest frequency to detect peaks.
%     bleachingLowpassFrequency - Lowest frequency to detect bleaching decay.
%     thresholdingFunction - One of @mad, @std.
%     thresholdFactor - Threshold cut-off.
%     triggeredWindow - Length of time to capture around each peak of spontaneous activity.
%     fitReference - Shift and scale reference to fit signal.
%     dffEpochs - Time epochs (s) to include for normalization.
%     f0Function - One of @movmean, @movmedian, @movmin.
%     f0Window - Length of the moving window to calculate f0.
%     f1Function - One of @movmean, @movmedian, @movmin, @movstd.
%     f1Window - Length of the moving window to calculate f1.
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
% results contain processed data, including df/f.
% 
% See examples
% See source code for detailed analysis steps and default parameters.
% Units for time and frequency are seconds and hertz respectively.

% 2019-02-01. Leonardo Molina.
% 2021-01-14. Last modified.
function results = FPA(time, signal, reference, configuration)
    results.warnings = {};
    if nargin < 4
        configuration = struct();
    end
    
    % Read input configuration. Use defaults for missing parameters.
    parameters.conditionEpochs = {'Data', [-Inf, Inf]};
    parameters.bleachingEpochs = [-Inf, Inf];
    parameters.artifactEpochs = [];
    parameters.resamplingFrequency = NaN;
    parameters.lowpassFrequency = 2;
    parameters.peaksLowpassFrequency = 0.5;
    parameters.bleachingLowpassFrequency = 0.1;
    parameters.fitReference = true;
    parameters.dffEpochs = [];
    parameters.f0Function = @movmedian;
    parameters.f0Window = Inf;
    parameters.f1Function = @movstd;
    parameters.f1Window = Inf;
    parameters.thresholdingFunction = @mad;
    parameters.thresholdFactor = 2.91;
    parameters.triggeredWindow = 10;
    parameters.f0 = [];
    parameters.f1 = [];
    parameters.plot = true;
    
    parameterNames = fieldnames(parameters);
    configurationNames = fieldnames(configuration);
    for i = 1:numel(configurationNames)
        name = configurationNames{i};
        if ismember(name, parameterNames)
            parameters.(name) = configuration.(name);
        else
            results.warnings{end + 1} = warn('[parsing] "%s" is not a valid parameter.', name);
        end
    end
    
    % Sampling frequency defaults.
    sourceFrequency = 1 / median(diff(time));
    if ~ismember('resamplingFrequency', configurationNames)
        parameters.resamplingFrequency = min(100, sourceFrequency);
    end
    
    % Plot: true | false | cell array of options.
    if islogical(parameters.plot)
        if parameters.plot
            parameters.plot = {'trace', 'power', 'stats', 'trigger'};
        else
            parameters.plot = {};
        end
    end
    
    % Settings for visualization.
    percentile = 0.99;
    grow = 0.50;
    
    % Make vectors of equal size.
    if isempty(reference)
        reference = NaN(size(signal));
        referenceProvided = false;
    elseif numel(reference) == 1 && isnan(reference)
        referenceProvided = false;
    else
        referenceProvided = true;
    end
    if referenceProvided
        k = isnan(time) | isnan(signal) | isnan(reference);
        time(k) = [];
        signal(k) = [];
        reference(k) = [];
    else
        k = isnan(time) | isnan(signal);
        time(k) = [];
        signal(k) = [];
    end
    
    % Resample to target frequency.
    if parameters.resamplingFrequency < sourceFrequency
        frequency = parameters.resamplingFrequency;
        % Express frequency as a ratio p/q.
        [p, q] = rat(frequency / sourceFrequency);
        % Resample: interpolate every p/q/f, upsample by p, filter, downsample by q.
        [signal, time2] = resample(signal, time, frequency, p, q);
        if referenceProvided
            reference = resample(reference, time, frequency, p, q);
        end
        time = time2;
    elseif parameters.resamplingFrequency ~= sourceFrequency
        frequency = sourceFrequency;
        results.warnings{end + 1} = warn('[resampling] Cannot resample to frequencies higher than the source frequency (%.2f Hz).', sourceFrequency);
    else
        frequency = sourceFrequency;
    end
    nSamples = numel(time);
    
    % Replace artifacts with straight lines.
    % Index of all points.
    allIds = colon(1, nSamples)';
    % Index of artifacts and non-artifacts.
    badId = time2id(time, parameters.artifactEpochs);
    goodId = setdiff(allIds, badId);
    % Interpolate.
    signal2 = signal;
    signal2(badId) = interp1(goodId, signal(goodId), badId);
    reference2 = reference;
    if referenceProvided
        reference2(badId) = interp1(goodId, reference(goodId), badId);
    end
    
    % Define clean epochs for fitting and peak detection.
    excludeWindow =  ceil(parameters.peaksLowpassFrequency * frequency);
    excludeIds = union(badId, [1:excludeWindow, numel(time) - excludeWindow + 1]');
    cleanIds = setdiff(time2id(time, cat(2, parameters.conditionEpochs{2:2:end})), excludeIds);
    
    % Model photo-bleaching with an exponential decay at given epochs.
    bleachingCorrectionId = time2id(time, parameters.bleachingEpochs);
    % Remove high-frequency oscillations to detect bleaching decay.
    if parameters.bleachingLowpassFrequency <= frequency / 2
        bleachingFilter = designfilt('lowpassiir', 'HalfPowerFrequency', parameters.bleachingLowpassFrequency, 'SampleRate', frequency, 'DesignMethod', 'butter', 'FilterOrder', 12);
        sLowpass = NaN(size(signal2));
        sLowpass(bleachingCorrectionId) = filtfilt(bleachingFilter, signal2(bleachingCorrectionId));
        rLowpass = NaN(size(reference2));
        if referenceProvided
            rLowpass(bleachingCorrectionId) = filtfilt(bleachingFilter, reference2(bleachingCorrectionId));
        end
    else
        rLowpass = reference2;
        sLowpass = signal2;
        results.warnings{end + 1} = warn('[bleaching-correction] Cannot lowpass to frequencies larger than half of the sampling frequency (%.2f Hz).', frequency / 2);
    end
    
    % Correct photobleaching on signal and reference.
    sFit = fit(time(bleachingCorrectionId), sLowpass(bleachingCorrectionId), fittype('exp1'));
    sBleaching = sFit(time);
    sCorrected = signal2 - sBleaching;
    if referenceProvided
        rFit = fit(time(bleachingCorrectionId), rLowpass(bleachingCorrectionId), fittype('exp1'));
        rBleaching = rFit(time);
        rCorrected = reference2 - rBleaching;
    else
        rCorrected = zeros(size(signal2));
    end
    
    % Fit reference to signal.
    if referenceProvided && parameters.fitReference
        r2sFit = fit(rCorrected(cleanIds), sCorrected(cleanIds), fittype('poly1'), 'Robust', 'on');
        rCorrected = r2sFit.p1 * rCorrected + r2sFit.p2;
    end
    
    % Correct for movement artifacts.
    f = sCorrected - rCorrected;
    
    % Low-pass filter.
    lowpassFilter = designfilt('lowpassiir', 'HalfPowerFrequency', parameters.lowpassFrequency, 'SampleRate', frequency, 'DesignMethod', 'butter', 'FilterOrder', 12);
    fLowpass = filtfilt(lowpassFilter, f);
    
    % Normalize.
    if isempty(parameters.f0) && isempty(parameters.f1)
        dffIds = time2id(time, parameters.dffEpochs);
        if isempty(dffIds)
            % For each data point, compute a baseline from its neighbors.
            n0 = nanmin(round(parameters.f0Window * frequency), numel(fLowpass));
            f0 = parameters.f0Function(fLowpass, n0);
            n1 = nanmin(round(parameters.f1Window * frequency), numel(fLowpass));
            f1 = parameters.f1Function(fLowpass, n1);
        else
            % For all points, compute a common baseline from the given epochs.
            f0 = parameters.f0Function(fLowpass(dffIds), numel(dffIds), 'Endpoints', 'discard');
            f1 = parameters.f1Function(fLowpass(dffIds), numel(dffIds), 'Endpoints', 'discard');
        end
    else
        f0 = parameters.f0;
        f1 = parameters.f1;
        if isempty(f0)
            f1 = f0;
        elseif isempty(f1)
            f0 = f1;
        end
    end
    dff = (fLowpass - f0) ./ f1;
    
    % Low-pass filter to detect peaks.
    if parameters.peaksLowpassFrequency <= frequency / 2
        peaksFilter = designfilt('lowpassiir', 'HalfPowerFrequency', parameters.peaksLowpassFrequency, 'SampleRate', frequency, 'DesignMethod', 'butter', 'FilterOrder', 12);
        peaksLowpass = filtfilt(peaksFilter, dff);
    else
        results.warnings{end + 1} = warn('[peak detection] Cannot lowpass to frequencies larger than half of the sampling frequency (%.2f Hz).', frequency / 2);
        peaksLowpass = dff;
    end
    
    % Get peak threshold.
    peakThreshold = mean(peaksLowpass(cleanIds)) + parameters.thresholdFactor * parameters.thresholdingFunction(peaksLowpass(cleanIds));
    valleyThreshold = -peakThreshold;
    
    state = warning('Query', 'signal:findpeaks:largeMinPeakHeight');
    warning('Off', 'signal:findpeaks:largeMinPeakHeight');
    if any(peaksLowpass >= peakThreshold)
        [~, peaksId] = findpeaks(+peaksLowpass, 'MinPeakHeight', peakThreshold);
        peaksId = intersect(peaksId, cleanIds);
    else
        peaksId = [];
    end
    if any(-peaksLowpass >= peakThreshold)
        [~, valleysId] = findpeaks(-peaksLowpass, 'MinPeakHeight', peakThreshold);
        valleysId = intersect(valleysId, cleanIds);
    else
        valleysId = [];
    end
    warning(state.state, 'signal:findpeaks:largeMinPeakHeight');
    
    % Split peaks/traces by conditions.
    % Number of samples in a triggered window.
    window = round(parameters.triggeredWindow * frequency);
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
    nEpochs = numel(parameters.conditionEpochs) / 2;
    epochBool = cell(1, nEpochs);
    for e = 1:nEpochs
        % Epoch index.
        ids = time2id(time, parameters.conditionEpochs{2 * e});
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
        ids = time2id(time, parameters.conditionEpochs{2 * e});
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
    
    signalColor = [0 0.4470 0.7410];
    referenceColor = [0.8500 0.3250 0.0980];
    alternateColor = [0 0.6470 0.9410];
    peaksLineColor = [0.4660 0.6740 0.1880];
    peaksMarkerColor = [1, 0, 0];
    dashColor = [0, 0, 0];
    results.figures = [];
    
    anyMatch = @(choices, pattern) any(~cellfun(@isempty, regexp(choices, pattern, 'start', 'once')));
    if anyMatch(parameters.plot, '\<trace\>')
        results.figures(end + 1) = figure('name', 'FPA: df/f');
        
        % Plot signal/reference resampled and bleaching model.
        ax.raw = subplot(5, 1, 1);
        ax.raw.XTick = [];
        hold(ax.raw, 'all');
        yy = [signal(:); reference(:); sBleaching(:)];
        ylims = limits(yy, percentile, grow);
        plotEpochs(parameters.conditionEpochs, xlims, ylims, cmap, true);
        plot(ax.raw, time, signal, 'Color', signalColor, 'DisplayName', 'Signal');
        if referenceProvided
            plot(ax.raw, time, reference, 'Color', referenceColor, 'DisplayName', 'Reference');
        end
        plot(ax.raw, time, sBleaching, 'Color', dashColor, 'LineStyle', '--', 'DisplayName', 'Photobleaching');
        ylim(ax.raw, ylims);
        title(ax.raw, 'Raw data');
        legend(ax.raw, 'show');
        
        % Plot signal/reference resampled and corrected for bleaching.
        ax.corrected = subplot(5, 1, 2);
        ax.corrected.XTick = [];
        hold(ax.corrected, 'all');
        yy = [sCorrected(:); rCorrected(:);];
        ylims = limits(yy, percentile, grow);
        plotEpochs(parameters.conditionEpochs, xlims, ylims, cmap, false);
        plot(ax.corrected, time, sCorrected, 'Color', signalColor, 'DisplayName', 'Signal');
        if referenceProvided
            plot(ax.corrected, time, rCorrected, 'Color', referenceColor, 'DisplayName', 'Reference');
        end
        ylim(ax.corrected, ylims);
        title(ax.corrected, 'Photobleaching correction');
        legend(ax.corrected, 'show');
        
        % Plot f and lowpass f.
        ax.f = subplot(5, 1, 3);
        ax.f.XTick = [];
        hold(ax.f, 'all');
        yy = [f(cleanIds); fLowpass(cleanIds)];
        ylims = limits(yy, percentile, grow);
        plotEpochs(parameters.conditionEpochs, xlims, ylims, cmap, false);
        plot(ax.f, time, f, 'Color', signalColor, 'DisplayName', 'f');
        plot(ax.f, time, fLowpass, 'Color', alternateColor, 'DisplayName', sprintf('f (<%.2fHz)', parameters.lowpassFrequency));
        ylim(ax.f, ylims);
        title(ax.f, 'Baseline correction');
        legend(ax.f, 'show');
        
        % Plot df/f.
        ax.filtered = subplot(5, 1, 4);
        ax.filtered.XTick = [];
        hold(ax.filtered, 'all');
        yy = [dff(cleanIds); peaksLowpass(cleanIds)];
        ylims = limits(yy, percentile, grow);
        epochs = parameters.conditionEpochs;
        epochs(1:2:end) = arrayfun(@(e) sprintf('area:%.2f', area(e)), 1:nEpochs, 'UniformOutput', false);
        plotEpochs(epochs, xlims, ylims, cmap, true);
        plot(ax.filtered, time, dff, 'Color', signalColor, 'DisplayName', 'df/f');
        plot(ax.filtered, time, peaksLowpass, 'Color', peaksLineColor, 'DisplayName', sprintf('df/f (<%.2fHz)', parameters.peaksLowpassFrequency));
        ylim(ax.filtered, ylims);
        title(ax.filtered, 'Normalization');
        legend(ax.filtered, 'show');

        % Plot df/f and peaks.
        ax.processed = subplot(5, 1, 5);
        hold(ax.processed, 'all');
        yy = peaksLowpass(cleanIds);
        ylims = limits(yy, percentile, grow);
        epochs = parameters.conditionEpochs;
        epochs(1:2:end) = arrayfun(@(e) sprintf('%i peaks / %i valleys', peaksCount(e), valleysCount(e)), 1:nEpochs, 'UniformOutput', false);
        plotEpochs(epochs, xlims, ylims, cmap, true);
        plot(ax.processed, time, peaksLowpass, 'Color', peaksLineColor, 'DisplayName', sprintf('df/f (<%.2fHz)', parameters.peaksLowpassFrequency));
        plot(ax.processed, time([1, end]), peakThreshold([1, 1]), 'Color', dashColor, 'LineStyle', '--', 'DisplayName', 'threshold');
        plot(ax.processed, time([1, end]), valleyThreshold([1, 1]), 'Color', dashColor, 'LineStyle', '--', 'HandleVisibility', 'off');
        plot(ax.processed, time(peaksId), peaksLowpass(peaksId), 'Color', peaksMarkerColor, 'LineStyle', 'none', 'Marker', 'o', 'HandleVisibility', 'off');
        plot(ax.processed, time(valleysId), peaksLowpass(valleysId), 'Color', peaksMarkerColor, 'LineStyle', 'none', 'Marker', 'o', 'HandleVisibility', 'off');
        ylim(ylims);
        title(ax.processed, 'Peak detection');
        legend(ax.processed, 'show');

        % Move axes together.
        linkaxes(findobj(gcf(), 'type', 'axes'), 'x');
        xlim(ax.raw, [time(1), time(end)]);

        xlabel('Time (s)');
        ylabel('df/f');
    end
    
    if anyMatch(parameters.plot, '\<power\>')
        % Plot power spectrum.
        results.figures(end + 1) = figure('name', 'FPA: Power spectrum');
        axs = cell(1, nEpochs);
        for e = 1:nEpochs
            axs{e} = subplot(nEpochs, 1, e);
            epochName = parameters.conditionEpochs{2 * e - 1};
            ids = time2id(time, parameters.conditionEpochs{2 * e});
            n = numel(ids);
            if n > 2
                d = dff(ids);
                halfN = floor(n / 2);
                f = fft(d);
                % Two-sided spectrum.
                p2 = abs(f / n);
                % Single-sided amplitude spectrum.
                p1 = p2(1:halfN + 1);
                p1(2:end-1) = 2 * p1(2:end-1);
                % Create frequency vector for range.
                fs = frequency * (0:halfN) / n;
                plot(fs, p1);
                ylim(limits(p1, percentile, grow));
            end
            title(sprintf('%s - Power spectrum', epochName));
        end
        ylabel('Power');
        xlabel('Frequency (Hz)');
        linkaxes(findobj(gcf(), 'type', 'axes'), 'x');
    end
    
    if anyMatch(parameters.plot, '\<stats\>')
        % Boxplot of dff.
        results.figures(end + 1) = figure('name', 'FPA: Boxplot');
        h = axes();
        % Not all epochs may be available.
        epochNames = parameters.conditionEpochs(1:2:end);
        groups = unique(epochGroups);
        boxplot(dff(epochIds), epochGroups, 'Labels', epochNames(groups));
        hold('all');
        area = zeros(1, nEpochs);
        ylims = ylim();
        for e = 1:nEpochs
            ids = time2id(time, parameters.conditionEpochs{2 * e});
            n = numel(ids);
            if n > 2
                area(e) = mean(dff(ids));
                text(e, ylims(2), epochStatLabels{e}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'Top');
            end
        end
        plot(h.XTick, area(groups), 'k.', 'DisplayName', 'Mean');
        ylabel('df/f');
        xtickangle(45);
        title('Stats on df/f traces for each condition');
    end
    
    if anyMatch(parameters.plot, '\<trigger\>')
        % Plot triggered average.
        results.figures(end + 1) = figure('name', 'FPA: Triggered average');
        ax.trigger = axes();
        hold(ax.trigger, 'all');
        if nGroups > 0
            timeTemplate = windowTemplate / frequency;
            for e = 1:nGroups
                group = uniqueGroups(e);
                triggeredDff = dff(triggeredId(peakGroups == group, :));
                triggeredDff = reshape(triggeredDff, numel(triggeredDff) / window, window);
                triggeredMean = mean(triggeredDff, 1);
                h1 = plot(timeTemplate, triggeredMean, 'HandleVisibility', 'off');
                epochName = parameters.conditionEpochs{2 * e - 1};
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
            text(ax.trigger, 0.5, 0.5, sprintf('No peaks above threshold %.2f (factor:%.2f)', peakThreshold, parameters.thresholdFactor), 'HorizontalAlignment', 'center');
        end
        title('Triggered average');
        legend('show');
        xlabel('Time (s)');
        ylabel('df/f');
        axis(ax.trigger, 'tight');
    end

    results.time = time;
    results.reference = reference2;
    results.signal = signal2;
    results.f = f;
    results.fLowpass = fLowpass;
    results.f0 = f0;
    results.f1 = f1;
    results.dff = dff;
    results.peaksLowpass = peaksLowpass;
    results.peaksId = peaksId;
    results.valleysId = valleysId;
    results.epochIds = epochIds;
    results.epochGroups = epochGroups;
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