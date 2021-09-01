% fpa = FPA(time, signal, reference, configuration);
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
%   -Plot 4, 5, 6:
%     -Peak triggers, epoch start-triggers, and epoch stop-triggers.
%   -Plot 7, 8, 9:
%     -Triggered average of peaks (if any), epoch start, and epoch stop.
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
% Fluorescence deflections are considered peaks when they exceed a threshold calculated as
% k * f2 + f3 and they are provided by the user as configuration.threshold = {k, f2, f3}
% Examples:
%   2.91 median absolute deviations from the median:
%     configuration.threshold = {2.91, @mad, @median}
%   2.91 median absolute deviations from 0:
%     configuration.threshold = {2.91, @mad, 0}
%   2.00 standard deviations from the mean:
%     configuration.threshold = {2.00, @std, @mean}
% 
% See examples and source code for detailed analysis steps and default parameters.
% Units for time and frequency are seconds and hertz respectively.
% 
% 2019-02-01. Leonardo Molina.
% 2021-09-01. Last modified.
classdef FPA < handle
    properties
        configuration
        
        warnings
        time
        frequency
        peakIds
        peakLabels
        epochIds
        epochLabels
        peakCounts
        valleyCounts
        
        signalRaw
        referenceRaw
        signalBaseline
        referenceBaseline
        signal
        reference
        f
        fSmooth
        f0
        f1
        dff
        area
        duration
        peakThreshold
    end
    
    properties (Access = private)
        nConditions
        epochNames
        cleanIds
        windowTemplate
        referenceProvided
        peaksSmooth
        uniquePeakIds
        uniqueValleyIds
        normalizedArea
        boxplotIds
        boxplotLabels
        
        % Settings for visualization.
        cmap = lines();
        zoomSettings = {0.99, 0.50};
        signalColor = [0, 0.4470, 0.7410];
        referenceColor = [0.8500, 0.3250, 0.0980];
        alternateColor = [0, 0.6470, 0.9410];
        peaksLineColor = [0.4660, 0.6740, 0.1880];
        peaksMarkerColor = [1, 0, 0];
        dashColor = [0, 0, 0];
    end
    
    methods
        function obj = FPA(time, signal, reference, parameters)
            obj.warnings = {};
            if nargin < 4
                parameters = struct();
            end

            % Read input configuration. Use defaults for missing parameters.
            defaults.conditionEpochs = {'Data', [-Inf, Inf]};
            defaults.artifactEpochs = [];
            defaults.resamplingFrequency = NaN;
            defaults.baselineEpochs = [-Inf, Inf];
            defaults.baselineLowpassFrequency = 0.1;
            defaults.airPLS = false;
            defaults.lowpassFrequency = 5;
            defaults.peaksLowpassFrequency = 0.5;
            defaults.fitReference = true;
            defaults.f0 = @median;
            defaults.f1 = @mad;
            defaults.threshold = {2.91, @mad, @median};
            defaults.triggeredWindow = 10;
            
            configuration = defaults;
            defaultNames = fieldnames(defaults);
            parametersNames = fieldnames(parameters);
            for i = 1:numel(parametersNames)
                name = parametersNames{i};
                if ismember(name, defaultNames)
                    configuration.(name) = parameters.(name);
                else
                    obj.warnings{end + 1} = warn('[parsing] "%s" is not a valid parameter.', name);
                end
            end

            % Sampling frequency defaults.
            sourceFrequency = 1 / median(diff(time));
            if ~ismember('resamplingFrequency', parametersNames)
                configuration.resamplingFrequency = min(100, sourceFrequency);
            end

            if islogical(configuration.airPLS)
                useAirPLS = configuration.airPLS;
                configuration.airPLS = [5e9, 2, 0.1, 0.5, 50];
            else
                useAirPLS = true;
            end
            configuration.airPLS = num2cell(configuration.airPLS);

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
            if configuration.resamplingFrequency < sourceFrequency
                frequency = configuration.resamplingFrequency;
                % Express frequency as a ratio p/q.
                [p, q] = rat(frequency / sourceFrequency);
                % Resample: interpolate every p/q/f, upsample by p, filter, downsample by q.
                [signal, time2] = resample(signal, time, frequency, p, q);
                if referenceProvided
                    reference = resample(reference, time, frequency, p, q);
                end
                time = time2;
            elseif configuration.resamplingFrequency ~= sourceFrequency
                frequency = sourceFrequency;
                obj.warnings{end + 1} = warn('[resampling] Cannot resample to frequencies higher than the source frequency (%.2f Hz).', sourceFrequency);
            else
                frequency = sourceFrequency;
            end

            % Setup.
            nSamples = numel(time);
            nConditions = numel(configuration.conditionEpochs) / 2;
            epochNames = configuration.conditionEpochs(1:2:end);

            % Replace artifacts with straight lines.
            % Index of all points.
            allIds = colon(1, nSamples)';
            % Index of artifacts and non-artifacts.
            badId = time2id(time, configuration.artifactEpochs);
            goodId = setdiff(allIds, badId);
            % Interpolate.
            signal2 = signal;
            signal2(badId) = interp1(goodId, signal(goodId), badId);
            reference2 = reference;
            if referenceProvided
                reference2(badId) = interp1(goodId, reference(goodId), badId);
            end
            
            % Define clean epochs for fitting and peak detection.
            excludeWindow =  ceil(0.5 * frequency / configuration.peaksLowpassFrequency);
            excludeIds = union(badId, [1:excludeWindow, numel(time) - excludeWindow + 1:nSamples]');
            cleanIds = setdiff(allIds, excludeIds);
            baselineCorrectionId = time2id(time, configuration.baselineEpochs);

            % Remove high-frequency oscillations to detect baseline (where indicated).
            signalSmooth = signal2;
            referenceSmooth = reference2;
            if configuration.baselineLowpassFrequency <= frequency / 2
                baselineFilter = designfilt('lowpassiir', 'HalfPowerFrequency', configuration.baselineLowpassFrequency, 'SampleRate', frequency, 'DesignMethod', 'butter', 'FilterOrder', 12);
                signalSmooth(baselineCorrectionId) = filtfilt(baselineFilter, signal2(baselineCorrectionId));
                if referenceProvided
                    referenceSmooth(baselineCorrectionId) = filtfilt(baselineFilter, reference2(baselineCorrectionId));
                end
            else
                obj.warnings{end + 1} = warn('[baseline-correction] Cannot lowpass to frequencies larger than half of the sampling frequency (%.2f Hz).', frequency / 2);
            end

            if useAirPLS
                % Model baseline with airPLS (everywhere).
                [~, signalBaseline] = airPLS(signalSmooth', configuration.airPLS{:});
                signalBaseline = signalBaseline';
                signalCorrected = signal2 - signalBaseline;
                if referenceProvided
                    [~, referenceBaseline] = airPLS(referenceSmooth', configuration.airPLS{:});
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
            if referenceProvided && configuration.fitReference
                r2sFit = fit(referenceCorrected(baselineCorrectionId), signalCorrected(baselineCorrectionId), fittype('poly1'), 'Robust', 'on');
                referenceCorrected = r2sFit.p1 * referenceCorrected + r2sFit.p2;
            end

            % Correct for movement artifacts.
            f = signalCorrected - referenceCorrected;

            % Low-pass filter.
            fFilter = designfilt('lowpassiir', 'HalfPowerFrequency', configuration.lowpassFrequency, 'SampleRate', frequency, 'DesignMethod', 'butter', 'FilterOrder', 12);
            fSmooth = f;
            fSmooth(cleanIds) = filtfilt(fFilter, f(cleanIds));

            % Normalize.
            f0 = normalize(configuration.f0, fSmooth, time);
            f1 = normalize(configuration.f1, fSmooth, time);
            dff = (fSmooth - f0) ./ f1;

            % Low-pass filter to detect peaks.
            peaksSmooth = dff;
            if configuration.peaksLowpassFrequency <= frequency / 2
                peaksFilter = designfilt('lowpassiir', 'HalfPowerFrequency', configuration.peaksLowpassFrequency, 'SampleRate', frequency, 'DesignMethod', 'butter', 'FilterOrder', 12);
                peaksSmooth(cleanIds) = filtfilt(peaksFilter, dff(cleanIds));
            else
                obj.warnings{end + 1} = warn('[peak detection] Cannot lowpass to frequencies larger than half of the sampling frequency (%.2f Hz).', frequency / 2);
            end

            % Get peak threshold.
            peakThreshold = threshold(configuration.threshold, peaksSmooth(cleanIds));

            state = warning('Query', 'signal:findpeaks:largeMinPeakHeight');
            warning('Off', 'signal:findpeaks:largeMinPeakHeight');
            if any(peaksSmooth >= peakThreshold)
                [~, uniquePeakIds] = findpeaks(+peaksSmooth, 'MinPeakHeight', peakThreshold);
                uniquePeakIds = intersect(uniquePeakIds, cleanIds);
            else
                uniquePeakIds = [];
            end
            if any(-peaksSmooth >= peakThreshold)
                [~, uniqueValleyIds] = findpeaks(-peaksSmooth, 'MinPeakHeight', peakThreshold);
                uniqueValleyIds = intersect(uniqueValleyIds, cleanIds);
            else
                uniqueValleyIds = [];
            end
            warning(state.state, 'signal:findpeaks:largeMinPeakHeight');
            [~, halfWindow] = forceOdd(configuration.triggeredWindow * frequency);
            % Index template to apply around each peak.
            obj.windowTemplate = -halfWindow:halfWindow;

            % Index epochs.
            % Start and stop vector indices for all provided epochs.
            epochIds = zeros(2, 0);
            % Numeric label corresponding to each epoch range.
            epochLabels = zeros(0, 1);
            % Vector index for each peak.
            peakIds = zeros(0, 1);
            % Numeric label corresponding to each peak.
            peakLabels = zeros(0, 1);
            % Misc indexing / labeling.
            boxplotIds = zeros(0, 1);
            boxplotLabels = zeros(0, 1);
            peakCounts = zeros(nConditions, 1);
            valleyCounts = zeros(nConditions, 1);
            area = zeros(nConditions, 1);
            duration = zeros(nConditions, 1);

            for c = 1:nConditions
                % Epoch indices:  3 4 5 ...
                [ids, bounds] = time2id(time, configuration.conditionEpochs{2 * c});
                nLimits = numel(bounds) / 2;
                epochIds = cat(2, epochIds, bounds);
                epochLabels = cat(1, epochLabels, repmat(c, nLimits, 1));

                % Triggered windows and condition labels (overlapping is possible and allowed).
                % Peaks in epoch: 3   5 ...
                epochPeakIds = intersect(ids, uniquePeakIds);
                nPeaks = numel(epochPeakIds);
                peakIds = cat(1, peakIds, epochPeakIds);
                peakLabels = cat(1, peakLabels, repmat(c, nPeaks, 1));

                boxplotIds = cat(1, boxplotIds, ids);
                boxplotLabels = cat(1, boxplotLabels, repmat(c, numel(ids), 1));

                peakCounts(c) = sum(ismember(ids, uniquePeakIds));
                valleyCounts(c) = sum(ismember(ids, uniqueValleyIds));

                duration(c) = numel(ids) / frequency;
                area(c) = trapz(dff(ids)) / frequency;
            end
            normalizedArea = area ./ duration;
            normalizedArea(duration == 0) = 0;
            
            obj.configuration = configuration;
            obj.time = time;
            obj.frequency = frequency;
            
            % Order depends on epoch definitions. Overlapping is possible and allowed.
            obj.peakIds = peakIds;
            obj.peakLabels = peakLabels;
            obj.epochIds = epochIds;
            obj.epochLabels = epochLabels;
            obj.peakCounts = peakCounts;
            obj.valleyCounts = valleyCounts;
            
            % Resampled only.
            obj.signalRaw = signal;
            obj.referenceRaw = reference;
            
            % Filtered, uncorrected.
            obj.signalBaseline = signalBaseline;
            obj.referenceBaseline = referenceBaseline;

            % Unfiltered, corrected.
            obj.signal = signalCorrected;
            obj.reference = referenceCorrected;

            % Unfiltered, motion corrected.
            obj.f = f;

            % Filtered, motion corrected.
            obj.fSmooth = fSmooth;

            % From fSmooth.
            obj.f0 = f0;
            obj.f1 = f1;
            obj.dff = dff;
            obj.area = area;
            obj.duration = duration;
            
            obj.nConditions = nConditions;
            obj.epochNames = epochNames;
            obj.cleanIds = cleanIds;
            obj.peakThreshold = peakThreshold;
            obj.peaksSmooth = peaksSmooth;
            obj.uniquePeakIds = uniquePeakIds;
            obj.uniqueValleyIds = uniqueValleyIds;
            obj.boxplotIds = boxplotIds;
            obj.boxplotLabels = boxplotLabels;
            obj.normalizedArea = normalizedArea;
            obj.referenceProvided = referenceProvided;
        end
        
        function plot(obj)
            obj.plotTrace();
            obj.plotPowerSpectrum();
            obj.plotStatistics();
            obj.plotTrigger();
            obj.plotTriggerAverage();
            obj.plotAUC();
        end
        
        function fig = plotTrace(obj)
            fig = figure('name', 'FPA: df/f');
            xlims = obj.time([1, end]);

            % Plot raw signal, reference, and baseline model.
            subplot(5, 1, 1);
            hold('all');
            yy = [obj.signalRaw(:); obj.referenceRaw(:); obj.signalBaseline(:)];
            ylims = limits(yy, obj.zoomSettings{:});
            plotEpochs(obj.configuration.conditionEpochs, xlims, ylims, obj.cmap, true);
            plot(obj.time, obj.signalRaw, 'Color', obj.signalColor, 'DisplayName', 'Signal');
            if obj.referenceProvided
                plot(obj.time, obj.referenceRaw, 'Color', obj.referenceColor, 'DisplayName', 'Reference');
                plot(obj.time, obj.referenceBaseline, 'Color', obj.dashColor, 'LineStyle', '--', 'DisplayName', 'Baseline', 'HandleVisibility', 'off');
            end
            plot(obj.time, obj.signalBaseline, 'Color', obj.dashColor, 'LineStyle', '--', 'DisplayName', 'Baseline');
            ylim(ylims);
            title('Raw data');
            legend('show');

            % Plot baseline corrected signal and reference.
            subplot(5, 1, 2);
            hold('all');
            yy = [obj.signal(:); obj.reference(:)];
            ylims = limits(yy, obj.zoomSettings{:});
            plotEpochs(obj.configuration.conditionEpochs, xlims, ylims, obj.cmap, false);
            plot(obj.time, obj.signal, 'Color', obj.signalColor, 'DisplayName', 'Signal');
            if obj.referenceProvided
                plot(obj.time, obj.reference, 'Color', obj.referenceColor, 'DisplayName', 'Reference');
            end
            ylim(ylims);
            title('Baseline correction');
            legend('show');

            % Plot motion correction (f and lowpass f).
            subplot(5, 1, 3);
            hold('all');
            yy = [obj.f(obj.cleanIds); obj.fSmooth(obj.cleanIds)];
            ylims = limits(yy, obj.zoomSettings{:});
            plotEpochs(obj.configuration.conditionEpochs, xlims, ylims, obj.cmap, false);
            plot(obj.time, obj.f, 'Color', obj.signalColor, 'DisplayName', 'f');
            plot(obj.time, obj.fSmooth, 'Color', obj.alternateColor, 'DisplayName', sprintf('f (<%.2fHz)', obj.configuration.lowpassFrequency));
            ylim(ylims);
            title('Motion correction');
            legend('show');

            % Plot normalization (e.g. df/f).
            subplot(5, 1, 4);
            hold('all');
            yy = [obj.dff(obj.cleanIds); obj.peaksSmooth(obj.cleanIds)];
            ylims = limits(yy, obj.zoomSettings{:});
            epochs = obj.configuration.conditionEpochs;
            epochs(1:2:end) = arrayfun(@(e) sprintf('area:%.2f', obj.area(e)), 1:obj.nConditions, 'UniformOutput', false);
            plotEpochs(epochs, xlims, ylims, obj.cmap, true);
            plot(obj.time, obj.dff, 'Color', obj.signalColor, 'DisplayName', 'df/f');
            plot(obj.time, obj.peaksSmooth, 'Color', obj.peaksLineColor, 'DisplayName', sprintf('df/f (<%.2fHz)', obj.configuration.peaksLowpassFrequency));
            ylim(ylims);
            title('Normalization');
            legend('show');

            % Plot peak detection.
            subplot(5, 1, 5);
            hold('all');
            yy = obj.peaksSmooth(obj.cleanIds);
            ylims = limits(yy, obj.zoomSettings{:});
            epochs = obj.configuration.conditionEpochs;
            epochs(1:2:end) = arrayfun(@(e) sprintf('%i peaks / %i valleys', obj.peakCounts(e), obj.valleyCounts(e)), 1:obj.nConditions, 'UniformOutput', false);
            plotEpochs(epochs, xlims, ylims, obj.cmap, true);
            plot(obj.time, obj.peaksSmooth, 'Color', obj.peaksLineColor, 'DisplayName', sprintf('df/f (<%.2fHz)', obj.configuration.peaksLowpassFrequency));
            plot(obj.time([1, end]), +obj.peakThreshold([1, 1]), 'Color', obj.dashColor, 'LineStyle', '--', 'DisplayName', sprintf('threshold:%.2f', obj.peakThreshold));
            plot(obj.time([1, end]), -obj.peakThreshold([1, 1]), 'Color', obj.dashColor, 'LineStyle', '--', 'HandleVisibility', 'off');
            plot(obj.time(obj.uniquePeakIds), obj.peaksSmooth(obj.uniquePeakIds), 'Color', obj.peaksMarkerColor, 'LineStyle', 'none', 'Marker', 'o', 'HandleVisibility', 'off');
            plot(obj.time(obj.uniqueValleyIds), obj.peaksSmooth(obj.uniqueValleyIds), 'Color', obj.peaksMarkerColor, 'LineStyle', 'none', 'Marker', 'o', 'HandleVisibility', 'off');
            ylim(ylims);
            title('Peak detection');
            legend('show');

            % Move axes together.
            axs = findobj(gcf(), 'type', 'axes');
            linkaxes(axs, 'x');
            xlim(axs(1), [obj.time(1), obj.time(end)]);
            set(axs(2:end), 'XTick', []);
            xlabel('Time (s)');
            ylabel('df/f');
        end
        
        function fig = plotPowerSpectrum(obj)
            % Plot power spectrum.
            fig = figure('name', 'FPA: Power spectrum');
            axs = cell(1, obj.nConditions);
            for c = 1:obj.nConditions
                axs{c} = subplot(obj.nConditions, 1, c);
                epochName = obj.configuration.conditionEpochs{2 * c - 1};
                ids = time2id(obj.time, obj.configuration.conditionEpochs{2 * c});
                n = numel(ids);
                if n > 2
                    d = obj.dff(ids);
                    halfN = floor(n / 2);
                    ff = fft(d);
                    % Two-sided spectrum.
                    p2 = abs(ff / n);
                    % Single-sided amplitude spectrum.
                    p1 = p2(1:halfN + 1);
                    p1(2:end-1) = 2 * p1(2:end-1);
                    % Create frequency vector for range.
                    fs = obj.frequency * (0:halfN) / n;
                    plot(fs, p1);
                    ylim(limits(p1, obj.zoomSettings{:}));
                end
                title(sprintf('%s - Power spectrum', epochName));
            end
            ylabel('Power');
            xlabel('Frequency (Hz)');
            linkaxes(findobj(gcf(), 'type', 'axes'), 'x');
        end

        function fig = plotStatistics(obj)
            % Boxplot of dff.
            fig = figure('name', 'FPA: Boxplot');
            % Not all epochs may be available.
            boxplotNames = obj.epochNames(unique(obj.boxplotLabels));
            boxplot(obj.dff(obj.boxplotIds), obj.boxplotLabels, 'Labels', boxplotNames);
            hold('all');
            ylims = ylim();
            for c = 1:obj.nConditions
                ids = time2id(obj.time, obj.configuration.conditionEpochs{2 * c});
                n = numel(ids);
                if n > 2
                    epochStatLabel = sprintf('\nmean:%.2f\nstd:%.2f', obj.normalizedArea(c), std(obj.dff(ids)));
                    text(c, ylims(2), epochStatLabel, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'Top');
                end
            end
            plot(obj.normalizedArea, 'k.', 'DisplayName', 'Mean');
            ylabel('df/f');
            xtickangle(45);
            title('Stats on df/f traces for each condition');
        end
            
        function figs = plotTrigger(obj)
            name = 'Peak-trigger heatmap';
            figs(1) = figure('name', name);
            plotTriggerHeatmap(obj.dff, obj.peakIds, obj.peakLabels, obj.windowTemplate, obj.frequency, obj.configuration.conditionEpochs(1:2:end));
            annotation('textbox', [0, 0.95, 1, 0.05], 'string', name, 'LineStyle', 'none');

            name = 'Start-triggered heatmap';
            figs(2) = figure('name', name);
            plotTriggerHeatmap(obj.dff, obj.epochIds(1:2:end)', obj.epochLabels, obj.windowTemplate, obj.frequency, obj.configuration.conditionEpochs(1:2:end));
            annotation('textbox', [0, 0.95, 1, 0.05], 'string', name, 'LineStyle', 'none');

            name = 'Stop-triggered heatmap';
            figs(3) = figure('name', name);
            plotTriggerHeatmap(obj.dff, obj.epochIds(2:2:end)', obj.epochLabels, obj.windowTemplate, obj.frequency, obj.configuration.conditionEpochs(1:2:end));
            annotation('textbox', [0, 0.95, 1, 0.05], 'string', name, 'LineStyle', 'none');
        end

        function figs = plotTriggerAverage(obj)
            figs(1) = figure('name', 'FPA: Peak-triggered average');
            plotTriggerAverage(obj.dff, obj.peakIds, obj.peakLabels, obj.windowTemplate, obj.frequency, obj.configuration.conditionEpochs(1:2:end), obj.cmap);
            ylabel('df/f');
            title('Peak-triggered average');

            figs(2) = figure('name', 'FPA: start-triggered average');
            plotTriggerAverage(obj.dff, obj.epochIds(1:2:end)', obj.epochLabels, obj.windowTemplate, obj.frequency, obj.configuration.conditionEpochs(1:2:end), obj.cmap);
            ylabel('df/f');
            title('Start-triggered average');

            figs(3) = figure('name', 'FPA: stop-triggered average');
            plotTriggerAverage(obj.dff, obj.epochIds(2:2:end)', obj.epochLabels, obj.windowTemplate, obj.frequency, obj.configuration.conditionEpochs(1:2:end), obj.cmap);
            ylabel('df/f');
            title('Stop-triggered average');
        end

        function fig = plotAUC(obj)
            % Plot normalized area under the curve.
            obj.nConditions = numel(obj.configuration.conditionEpochs) / 2;
            fig = figure('name', 'FPA: Normalized area under the curve');
            ax.auc = axes();
            bar(1:obj.nConditions, obj.normalizedArea);
            set(ax.auc, 'XTickLabel', obj.epochNames);
            xtickangle(45);
            title('dff/f - normalized AUC');
        end
        
        function export(obj, prefix)
            prefix = regexprep(prefix, '/$', '');
            [folder, basename] = fileparts(prefix);
            
            % Save data for post-processing.

            % Time vs dff.
            % Rows represent increasing values of time with corresponding dff values.
            output = fullfile(folder, sprintf('%s - dff.csv', basename));
            fid = fopen(output, 'w');
            fprintf(fid, '# time, dff\n');
            fprintf(fid, '%.4f, %.4f\n', [obj.time, obj.dff]');
            fclose(fid);
            
            % AUC.
            output = fullfile(folder, sprintf('%s - AUC.csv', basename));
            fid = fopen(output, 'w');
            fprintf(fid, '# conditionId, conditionName, area, duration\n');
            data = [obj.epochNames(:), num2cell([colon(1, obj.nConditions)', obj.area, obj.duration])];
            data = data(:, [2, 1, 3, 4])';
            fprintf(fid, '%i, %s, %.4f, %d\n', data{:});
            fclose(fid);

            % All peak-triggered windows and their average with corresponding epoch label.
            output1 = fullfile(folder, sprintf('%s - peak-triggered.csv', basename));
            output2 = fullfile(folder, sprintf('%s - peak-triggered averaged.csv', basename));
            obj.saveEventTrigger(output1, output2, 'peak', obj.peakIds, obj.peakLabels);

            % Same as above for start-triggered windows.
            output1 = fullfile(folder, sprintf('%s - start-triggered.csv', basename));
            output2 = fullfile(folder, sprintf('%s - start-triggered averaged.csv', basename));
            startIds = obj.epochIds(1:2:end);
            obj.saveEventTrigger(output1, output2, 'condition', startIds, obj.epochLabels, obj.epochNames(obj.epochLabels));

            % Same as above for stop-triggered windows.
            output1 = fullfile(folder, sprintf('%s - stop-triggered.csv', basename));
            output2 = fullfile(folder, sprintf('%s - stop-triggered averaged.csv', basename));
            stopIds = obj.epochIds(2:2:end);
            obj.saveEventTrigger(output1, output2, 'condition', stopIds, obj.epochLabels, obj.epochNames(obj.epochLabels));
        end
    end
    
    methods (Access = private)
        function saveEventTrigger(obj, output1, output2, prefix, triggerIds, triggerLabels, triggerNames)
            [triggeredWindow, halfWindow] = forceOdd(obj.configuration.triggeredWindow * obj.frequency);
            window = -halfWindow:halfWindow;
            traceTime = window / obj.frequency;
            timeHeader = strjoin(arrayfun(@(x) sprintf('%.4f', x), traceTime, 'UniformOutput', false), ', ');
            
            % Filter out out-of-range traces.
            nSamples = numel(obj.dff);
            k = triggerIds > halfWindow & triggerIds + halfWindow < nSamples;
            triggerLabels = triggerLabels(k);
            triggerIds = triggerIds(k);
            
            triggerIds = triggerIds(:);
            triggerLabels = triggerLabels(:);
            if nargin == 7
                triggerNames = triggerNames(:);
                saveConditionName = true;
            else
                saveConditionName = false;
            end
            if numel(triggerIds) >= 1
                triggerData = obj.dff;
                triggerData = triggerData(window + triggerIds);
                triggerData = reshape(triggerData, [], numel(window));
                triggerTime = obj.time(triggerIds);

                % All triggered windows with corresponding epoch label.
                % Order depends on epoch definitions. Overlapping is possible.
                % Rows represent a single trigger:
                % First column is the condition label of the trigger followed by the trace around each trigger, with each trigger at the center column (n / 2 + 1) labeled with "zero".
                fid = fopen(output1, 'w');
                if saveConditionName
                    fprintf(fid, ['# time, %1$sId, %1$sName, ', timeHeader, '\n'], prefix);
                    format = ['%.4f, %i, %s', repmat(', %.4f', 1, triggeredWindow), '\n'];
                    data = [triggerNames, num2cell([triggerTime, triggerLabels, triggerData])];
                    data = data(:, [2 3 1 4:end])';
                else
                    fprintf(fid, ['# time, %1$sId, ', timeHeader, '\n'], prefix);
                    format = ['%.4f, %i', repmat(', %.4f', 1, triggeredWindow), '\n'];
                    data = num2cell([triggerTime, triggerLabels, triggerData]);
                    data = data(:, :)';
                end
                fprintf(fid, format, data{:});
                fclose(fid);

                % Average of the above.
                uLabels = unique(triggerLabels);
                averages = zeros(0, size(triggerData, 2));
                for u = 1:numel(uLabels)
                    label = uLabels(u);
                    uData = triggerData(triggerLabels == label, :);
                    uData = reshape(uData, [], size(triggerData, 2));
                    averages = cat(1, averages, mean(uData, 1));
                end
                
                fid = fopen(output2, 'w');
                if saveConditionName
                    fprintf(fid, ['# %1$sId, %1$sName, ', timeHeader, '\n'], prefix);
                    format = ['%i, %s', repmat(', %.4f', 1, triggeredWindow), '\n'];
                    uTriggerNames = unique(triggerNames, 'stable');
                    data = [uTriggerNames, num2cell([uLabels, averages])];
                    data = data(:, [2, 1, 3:end])';
                    fprintf(fid, format, data{:});
                else
                    fprintf(fid, ['# %1$sId, ', timeHeader, '\n'], prefix);
                    format = ['%i', repmat(', %.4f', 1, triggeredWindow), '\n'];
                    data = num2cell([uLabels, averages]);
                    data = data(:, :)';
                    fprintf(fid, format, data{:});
                end
                fclose(fid);
            end
        end
    end
end

function output = normalize(parameters, f, time)
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

function value = threshold(parameters, data)
    % {value1, @mad, @median}
    % {value1, @mad}
    % {value1, @mad, value2}
    % {value1}
    % value
    if iscell(parameters)
        n = numel(parameters);
        if n >= 1
            k = parameters{1};
        else
            k = 2.91;
        end
        if n >= 2
            f2 = parameters{2};
        else
            f2 = @mad;
        end
        if n >= 3
            f3 = parameters{3};
        else
            f3 = @median;
        end
    else
        k = parameters;
        f2 = @mad;
        f3 = @median;
    end
    if isa(f3, 'function_handle')
        value = k * f2(data) + f3(data);
    else
        value = k * f2(data) + f3;
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
    x = x(:);
    ylims = [prctile(x, 100 * (1 - percentile)), prctile(x, 100 * percentile)];
    delta = diff(ylims) * grow;
    ylims = [ylims(1) - delta, ylims(2) + delta];
end

function [odd, half] = forceOdd(fractional)
    % Number of samples in a triggered window (whole length, left to right).
    odd = fractional;
    % Force odd count.
    odd = round(odd) + (mod(round(odd), 2) == 0);
    half = (odd - 1) / 2;
end

function plotTriggerHeatmap(data, ids, labels, window, frequency, names)
    % Filter out out-of-range traces.
    nSamples = numel(data);
    triggeredWindow = numel(window);
    halfWindow = (triggeredWindow - 1) / 2;
    k = ids > halfWindow & ids + halfWindow < nSamples;
    labels = labels(k);
    ids = ids(k);
    time = window / frequency;
    
    nConditions = numel(names);
    [ni, nj] = squaredFactors(nConditions);
    
    if numel(ids) > 0
        allWindowIds = ids + window;
        allTriggeredData = data(allWindowIds);
        clims = limits(allTriggeredData, 0.99, 0);
    end
    
    axs = cell(nConditions, 1);
    for c = 1:nConditions
        axs{c} = subplot(ni, nj, c);
        triggerIds = ids(labels == c);
        nTriggers = numel(triggerIds);
        if nTriggers > 0
            windowIds = triggerIds + window;
            triggeredData = data(windowIds);
            % Make sure matrix is nr x nc, particularly for 1 x nc.
            triggeredData = reshape(triggeredData, size(windowIds));
            nTriggers = size(triggeredData, 1);
            imagesc('xData', time, 'yData', 1:nTriggers, 'cData', triggeredData, clims);
            yticks = get(gca(), 'YTick');
            yticks = yticks(round(yticks) == yticks);
            set(gca(), 'YTick', yticks);
        else
            text(0.5, 0.5, 'No triggers', 'HorizontalAlignment', 'center');
        end
        title(names{c});
        axis('tight');
    end
    xlabel(axs{ni}, 'Time (s)');
    ylabel(axs{ni}, 'Trigger id');
    set(get(colorbar(), 'title'), 'string', 'df/f');
end

function plotTriggerAverage(data, ids, labels, window, frequency, names, colors)
    % Filter out out-of-range traces.
    nSamples = numel(data);
    triggeredWindow = numel(window);
    halfWindow = (triggeredWindow - 1) / 2;
    k = ids > halfWindow & ids + halfWindow < nSamples;
    labels = labels(k);
    ids = ids(k);
    time = window / frequency;
    
    if numel(ids) > 0
        hold('all');
        nConditions = numel(names);
        for c = 1:nConditions
            triggerIds = ids(labels == c);
            nTriggers = numel(triggerIds);
            if nTriggers > 0
                windowIds = triggerIds + window;
                triggeredData = data(windowIds);
                % Make sure matrix is nr x nc, particularly for 1 x nc.
                triggeredData = reshape(triggeredData, size(windowIds));
                av = mean(triggeredData, 1);
                nTriggers = size(triggeredData, 1);
                % Plot.
                plot(time, av, 'Color', colors(c, :), 'HandleVisibility', 'off');
                sem = std(triggeredData, [], 1) / sqrt(size(triggeredData, 1));
                sem0 = sem(ceil(size(triggeredData, 2) / 2));
                label = sprintf('%s (SEM=%.4f, n = %i)', names{c}, sem0, nTriggers);
                vertices = [time; av + sem / 2];
                vertices = cat(2, vertices, [fliplr(time); fliplr(av - sem / 2)])';
                faces = 1:2 * triggeredWindow;
                patch('Faces', faces, 'Vertices', vertices, 'FaceColor', colors(c, :), 'EdgeColor', 'none', 'FaceAlpha', 0.10, 'DisplayName', label);
            end
        end
    else
        text(0.5, 0.5, 'No triggers', 'HorizontalAlignment', 'center');
    end
    legend('show');
    xlabel('Time (s)');
    axis('tight');
end