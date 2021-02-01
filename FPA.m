% results = FPA(time, signal, reference, configuration)
% 
% Remove baseline from data, correct from motion artifacts; normalize, filter,
% and detect peaks of spontaneous activity in user defined epochs.
% 
% Time, signal, and reference are column vectors.
% 
% Processing steps:
%   -Resample signal and reference to target frequency.
%   -Replace artifacts with linear interpolation in flagged regions.
% 	-Baseline correction modeled as an exponential decay of the low-pass
%    filtered data or using airPLS.
% 	-Correct for motion artifacts by subtracting reference to signal, after a polynomial fit.
%   -Remove fast oscillations with a low-pass filter.
% 	-Normalize data as df/f or z-score according to settings.
% 	-Find peaks of spontaneous activity in low-pass filtered data.
%   -Plot 1:
%     -Raw signal, reference, and baseline model.
%     -Baseline corrected signal and reference.
%     -Motion correction.
%     -Normalization.
%     -Peak detection.
%   -Plot 2:
%     -Power spectrum for each epoch.
%   -Plot 3:
%     -Boxplot.
%   -Plot 4:
%     -Triggered average of spontaneous activity (if any peaks are found).
% 
% configuration is a struct with the following fields (defaults are used for missing fields):
%     conditionEpochs - Epochs for different conditions: {'epoch1', [start1, end1, start2, end2, ...], 'epoch2', ...}
%     baselineEpochs - Time epochs (s) to include for baseline correction.
%     baselineLowpassFrequency - Frequency representative of baseline.
%     airPLS - Baseline correction for all data using airPLS (true, false, or airPLS inputs).
%     artifactEpochs - Time epochs (s) to remove.
%     resamplingFrequency - Resampling frequency (Hz).
%     lowpassFrequency - Lowest frequency permitted in normalized signal.
%     peaksLowpassFrequency - Lowest frequency to detect peaks.
%     thresholdingFunction - @mad, @std, ...
%     thresholdFactor - Threshold cut-off.
%     triggeredWindow - Length of time to capture around each peak of spontaneous activity.
%     fitReference - Shift and scale reference to fit signal.
% 
% Normalization is calculated as (f - f0) / f1 where f0 and f1 can be data provided by the
% user or calculated using given functions:
% 
%     Normalization from given functions:
%         f0 and f1 are common to all datapoints and calculated from all data:
%             df/f:
%                 configuration.f0 = @mean;
%                 configuration.f1 = @mean;
%             z-score:
%                 configuration.f0 = @mean;
%                 configuration.f1 = @std;
%             z-score - alternative 1 (default):
%                 configuration.f0 = @median;
%                 configuration.f1 = @mad;
%             z-score - alternative 2:
%                 configuration.f0 = @median;
%                 configuration.f1 = @std;
% 
%         f0 and f1 are common to all data points and calculated at given epochs:
%             df/f:
%                 epochs = [0, 100, 500, 550, 1000, Inf]
%                 configuration.f0 = {@mean, epochs};
%                 configuration.f1 = {@mean, epochs};
% 
%         f0 and f1 are calculated for each data point based on a moving window:
%             df/f:
%                 window = 60;
%                 configuration.f0 = {@movmean, window};
%                 configuration.f1 = {@movmean, window};
%             (further combinations possible with @movmean, @movmedian, @movstd, @movmad, @mov...).
% 
%     Normalization from given data:
%         f0 = ones(size(time));
%         f1 = ones(size(time)) * 10;
%         configuration.f0 = f0;
%         configuration.f1 = f1;
% 
% results contain processed data.
% 
% See examples
% See source code for detailed analysis steps and default parameters.
% Units for time and frequency are seconds and hertz respectively.
% 
% 2019-02-01. Leonardo Molina.
% 2021-02-01. Last modified.
function results = FPA(time, signal, reference, configuration)
    results.warnings = {};
    if nargin < 4
        configuration = struct();
    end
    
    % Read input configuration. Use defaults for missing parameters.
    parameters.conditionEpochs = {'Data', [-Inf, Inf]};
    parameters.artifactEpochs = [];
    parameters.resamplingFrequency = NaN;
    parameters.baselineEpochs = [-Inf, Inf];
    parameters.baselineLowpassFrequency = 0.1;
    parameters.airPLS = false;
    parameters.lowpassFrequency = 2;
    parameters.peaksLowpassFrequency = 0.5;
    parameters.fitReference = true;
    parameters.f0 = @median;
    parameters.f1 = @mad;
    parameters.thresholdingFunction = @mad;
    parameters.thresholdFactor = 2.91;
    parameters.triggeredWindow = 10;
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
    
    % Settings for visualization.
    percentile = 0.99;
    grow = 0.50;
    
    % Sampling frequency defaults.
    sourceFrequency = 1 / median(diff(time));
    if ~ismember('resamplingFrequency', configurationNames)
        parameters.resamplingFrequency = min(100, sourceFrequency);
    end
    
    % Plot: true | false | cell array of choices.
    if islogical(parameters.plot)
        if parameters.plot
            parameters.plot = {'trace', 'power', 'stats', 'trigger'};
        else
            parameters.plot = {};
        end
    end
    
    if islogical(parameters.airPLS)
        useAirPLS = parameters.airPLS;
        parameters.airPLS = [5e9, 2, 0.1, 0.5, 50];
    else
        useAirPLS = true;
    end
    parameters.airPLS = num2cell(parameters.airPLS);
    
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
    cleanIds = setdiff(allIds, excludeIds);
    
    baselineCorrectionId = time2id(time, parameters.baselineEpochs);
    
    % Remove high-frequency oscillations to detect baseline (where indicated).
    signalSmooth = signal2;
    referenceSmooth = reference2;
    if parameters.baselineLowpassFrequency <= frequency / 2
        baselineFilter = designfilt('lowpassiir', 'HalfPowerFrequency', parameters.baselineLowpassFrequency, 'SampleRate', frequency, 'DesignMethod', 'butter', 'FilterOrder', 12);
        signalSmooth(baselineCorrectionId) = filtfilt(baselineFilter, signal2(baselineCorrectionId));
        if referenceProvided
            referenceSmooth(baselineCorrectionId) = filtfilt(baselineFilter, reference2(baselineCorrectionId));
        end
    else
        results.warnings{end + 1} = warn('[baseline-correction] Cannot lowpass to frequencies larger than half of the sampling frequency (%.2f Hz).', frequency / 2);
    end
    
    if useAirPLS
        % Model baseline with airPLS (everywhere).
        [~, signalBaseline] = airPLS(signalSmooth', parameters.airPLS{:});
        signalBaseline = signalBaseline';
        signalCorrected = signal2 - signalBaseline;
        if referenceProvided
            [~, referenceBaseline] = airPLS(referenceSmooth', parameters.airPLS{:});
            referenceBaseline = referenceBaseline';
            referenceCorrected = reference2 - referenceBaseline;
        else
            referenceBaseline = zeros(size(signalCorrected));
            referenceCorrected = zeros(size(signalCorrected));
        end
    else
        % Model baseline with an exponential decay at given epochs (where indicated).
        signalFit = fit(time(baselineCorrectionId), signalSmooth(baselineCorrectionId), fittype('exp1'));
        signalBaseline = signalFit(time);
        signalCorrected = signal2 - signalBaseline;
        if referenceProvided
            referenceFit = fit(time(baselineCorrectionId), referenceSmooth(baselineCorrectionId), fittype('exp1'));
            referenceBaseline = referenceFit(time);
            referenceCorrected = reference2 - referenceBaseline;
        else
            referenceBaseline = zeros(size(signalCorrected));
            referenceCorrected = zeros(size(signalCorrected));
        end
    end
    
    % Fit reference to signal (where indicated).
    if referenceProvided && parameters.fitReference
        r2sFit = fit(referenceCorrected(baselineCorrectionId), signalCorrected(baselineCorrectionId), fittype('poly1'), 'Robust', 'on');
        referenceCorrected = r2sFit.p1 * referenceCorrected + r2sFit.p2;
    end
    
    % Correct for movement artifacts.
    f = signalCorrected - referenceCorrected;
    
    % Low-pass filter.
    fFilter = designfilt('lowpassiir', 'HalfPowerFrequency', parameters.lowpassFrequency, 'SampleRate', frequency, 'DesignMethod', 'butter', 'FilterOrder', 12);
    fSmooth = f;
    fSmooth(cleanIds) = filtfilt(fFilter, f(cleanIds));
    
    % Normalize.
    f0 = parseNormalization(parameters.f0, fSmooth, time);
    f1 = parseNormalization(parameters.f1, fSmooth, time);
    dff = (fSmooth - f0) ./ f1;
    
    % Low-pass filter to detect peaks.
    peaksSmooth = dff;
    if parameters.peaksLowpassFrequency <= frequency / 2
        peaksFilter = designfilt('lowpassiir', 'HalfPowerFrequency', parameters.peaksLowpassFrequency, 'SampleRate', frequency, 'DesignMethod', 'butter', 'FilterOrder', 12);
        peaksSmooth(cleanIds) = filtfilt(peaksFilter, dff(cleanIds));
    else
        results.warnings{end + 1} = warn('[peak detection] Cannot lowpass to frequencies larger than half of the sampling frequency (%.2f Hz).', frequency / 2);
    end
    
    % Get peak threshold.
    peakThreshold = mean(peaksSmooth(cleanIds)) + parameters.thresholdFactor * parameters.thresholdingFunction(peaksSmooth(cleanIds));
    valleyThreshold = -peakThreshold;
    
    state = warning('Query', 'signal:findpeaks:largeMinPeakHeight');
    warning('Off', 'signal:findpeaks:largeMinPeakHeight');
    if any(peaksSmooth >= peakThreshold)
        [~, peaksId] = findpeaks(+peaksSmooth, 'MinPeakHeight', peakThreshold);
        peaksId = intersect(peaksId, cleanIds);
    else
        peaksId = [];
    end
    if any(-peaksSmooth >= peakThreshold)
        [~, valleysId] = findpeaks(-peaksSmooth, 'MinPeakHeight', peakThreshold);
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
    for e = 1:nEpochs
        % Epoch indices:     3 4 5 ... nSamples
        ids = time2id(time, parameters.conditionEpochs{2 * e});
        % Epoch flagged: 0 0 1 1 1 ... nSamples
        epochBool = ismember(allIds, ids);
        % Flag peaks within epoch: [0 0 1 1 1][peaksId]
        k = epochBool(peaksId);
        % Warn about overlaps.
        overlaps = peakGroups(k);
        overlaps = unique(overlaps(overlaps > 0));
        epochName = parameters.conditionEpochs{2 * e - 1};
        for o = overlaps(:)'
            overlapName = parameters.conditionEpochs{2 * o - 1};
            results.warnings{end + 1} = warn('[peaks] "%s" overlaps with "%s".', epochName, overlapName);
        end
        % Assign group to such peaks. When epochs overlap, the later overrides the group .
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
        
        % Plot raw signal, reference, and baseline model.
        ax.raw = subplot(5, 1, 1);
        ax.raw.XTick = [];
        hold(ax.raw, 'all');
        yy = [signal(:); reference(:); signalBaseline(:)];
        ylims = limits(yy, percentile, grow);
        plotEpochs(parameters.conditionEpochs, xlims, ylims, cmap, true);
        plot(ax.raw, time, signal, 'Color', signalColor, 'DisplayName', 'Signal');
        if referenceProvided
            plot(ax.raw, time, reference, 'Color', referenceColor, 'DisplayName', 'Reference');
        end
        plot(ax.raw, time, signalBaseline, 'Color', dashColor, 'LineStyle', '--', 'DisplayName', 'Baseline');
        ylim(ax.raw, ylims);
        title(ax.raw, 'Raw data');
        legend(ax.raw, 'show');
        
        % Plot baseline corrected signal and reference.
        ax.corrected = subplot(5, 1, 2);
        ax.corrected.XTick = [];
        hold(ax.corrected, 'all');
        yy = [signalCorrected(:); referenceCorrected(:);];
        ylims = limits(yy, percentile, grow);
        plotEpochs(parameters.conditionEpochs, xlims, ylims, cmap, false);
        plot(ax.corrected, time, signalCorrected, 'Color', signalColor, 'DisplayName', 'Signal');
        if referenceProvided
            plot(ax.corrected, time, referenceCorrected, 'Color', referenceColor, 'DisplayName', 'Reference');
        end
        ylim(ax.corrected, ylims);
        title(ax.corrected, 'Baseline correction');
        legend(ax.corrected, 'show');
        
        % Plot motion correction (f and lowpass f).
        ax.f = subplot(5, 1, 3);
        ax.f.XTick = [];
        hold(ax.f, 'all');
        yy = [f(cleanIds); fSmooth(cleanIds)];
        ylims = limits(yy, percentile, grow);
        plotEpochs(parameters.conditionEpochs, xlims, ylims, cmap, false);
        plot(ax.f, time, f, 'Color', signalColor, 'DisplayName', 'f');
        plot(ax.f, time, fSmooth, 'Color', alternateColor, 'DisplayName', sprintf('f (<%.2fHz)', parameters.lowpassFrequency));
        ylim(ax.f, ylims);
        title(ax.f, 'Motion correction');
        legend(ax.f, 'show');
        
        % Plot normalization (e.g. df/f).
        ax.filtered = subplot(5, 1, 4);
        ax.filtered.XTick = [];
        hold(ax.filtered, 'all');
        yy = [dff(cleanIds); peaksSmooth(cleanIds)];
        ylims = limits(yy, percentile, grow);
        epochs = parameters.conditionEpochs;
        epochs(1:2:end) = arrayfun(@(e) sprintf('area:%.2f', area(e)), 1:nEpochs, 'UniformOutput', false);
        plotEpochs(epochs, xlims, ylims, cmap, true);
        plot(ax.filtered, time, dff, 'Color', signalColor, 'DisplayName', 'df/f');
        plot(ax.filtered, time, peaksSmooth, 'Color', peaksLineColor, 'DisplayName', sprintf('df/f (<%.2fHz)', parameters.peaksLowpassFrequency));
        ylim(ax.filtered, ylims);
        title(ax.filtered, 'Normalization');
        legend(ax.filtered, 'show');

        % Plot peak detection.
        ax.processed = subplot(5, 1, 5);
        hold(ax.processed, 'all');
        yy = peaksSmooth(cleanIds);
        ylims = limits(yy, percentile, grow);
        epochs = parameters.conditionEpochs;
        epochs(1:2:end) = arrayfun(@(e) sprintf('%i peaks / %i valleys', peaksCount(e), valleysCount(e)), 1:nEpochs, 'UniformOutput', false);
        plotEpochs(epochs, xlims, ylims, cmap, true);
        plot(ax.processed, time, peaksSmooth, 'Color', peaksLineColor, 'DisplayName', sprintf('df/f (<%.2fHz)', parameters.peaksLowpassFrequency));
        plot(ax.processed, time([1, end]), peakThreshold([1, 1]), 'Color', dashColor, 'LineStyle', '--', 'DisplayName', 'threshold');
        plot(ax.processed, time([1, end]), valleyThreshold([1, 1]), 'Color', dashColor, 'LineStyle', '--', 'HandleVisibility', 'off');
        plot(ax.processed, time(peaksId), peaksSmooth(peaksId), 'Color', peaksMarkerColor, 'LineStyle', 'none', 'Marker', 'o', 'HandleVisibility', 'off');
        plot(ax.processed, time(valleysId), peaksSmooth(valleysId), 'Color', peaksMarkerColor, 'LineStyle', 'none', 'Marker', 'o', 'HandleVisibility', 'off');
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
            for g = 1:nGroups
                e = uniqueGroups(g);
                triggeredDff = dff(triggeredId(peakGroups == e, :));
                triggeredDff = reshape(triggeredDff, numel(triggeredDff) / window, window);
                triggeredMean = mean(triggeredDff, 1);
                plot(timeTemplate, triggeredMean, 'Color', cmap(e, :), 'HandleVisibility', 'off');
                epochName = parameters.conditionEpochs{2 * e - 1};
                triggeredSem = std(triggeredDff, [], 1) / sqrt(size(triggeredDff, 1));
                semAtZero = triggeredSem(ceil(window / 2));
                nPeaks = sum(peakGroups == e);
                label = sprintf('%s (SEM=%.4f, n = %i)', epochName, semAtZero, nPeaks);
                vertices = [timeTemplate; triggeredMean + triggeredSem / 2];
                vertices = cat(2, vertices, [fliplr(timeTemplate); fliplr(triggeredMean - triggeredSem / 2)])';
                faces = 1:2 * window;
                patch('Faces', faces, 'Vertices', vertices, 'FaceColor', cmap(e, :), 'EdgeColor', 'none', 'FaceAlpha', 0.10, 'DisplayName', label);
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

    % Filtered, uncorrected.
    results.signalBaseline = signalBaseline;
    results.referenceBaseline = referenceBaseline;
    
    % Unfiltered, corrected.
    results.signal = signalCorrected;
    results.reference = referenceCorrected;
    
    % Unfiltered, motion corrected.
    results.f = f;
    
    % Filtered, motion corrected.
    results.fSmooth = fSmooth;
    
    % From fSmooth.
    results.f0 = f0;
    results.f1 = f1;
    results.dff = dff;
    results.area = area;
    
    % fSmooth filtered again.
    results.peaks = peaksSmooth;
    results.peaksId = peaksId;
    results.valleysId = valleysId;
end

function output = parseNormalization(parameters, f, time)
    if iscell(parameters)
        fcn = parameters{1};
        if numel(parameters) == 1
            parameters{2} = [-Inf, Inf];
        end
        if isscalar(parameters{2})
            % Produce a vector from moving window.
            if numel(parameters) <= 2
                options = {'EndPoints', 'shrink'};
            else
                options = parameters(3:end);
            end
            frequency = 1 / median(diff(time));
            nSamples = numel(time);
            window = parameters{2};
            window = min(round(window * frequency), nSamples);
            output = fcn(f, window, options{:});
        else
            % Produce a value from all data (or epochs).
            epochs = parameters{2};
            ids = time2id(time, epochs);
            output = fcn(f(ids));
        end
    elseif isa(parameters, 'function_handle')
        % Produce a value from all data (or epochs).
        fcn = parameters;
        epochs = [-Inf, Inf];
        ids = time2id(time, epochs);
        output = fcn(f(ids));
    else
        output = parameters;
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