% fpa = FPA(time, signal, reference, configuration);
% 
% Remove baseline from data, correct from motion artifacts; normalize, filter,
% and detect peaks of spontaneous activity in user defined epochs.
% 
% Time, signal, and reference are column vectors.
% 
% Processing steps:
%   -Resample signal and reference to a given frequency.
%   -Replace artifacts with linear interpolation in flagged regions.
% 	-Baseline correction modeled as an exponential decay of the low-pass filtered data (optionally using airPLS).
% 	-Correct for motion artifacts by subtracting reference to signal, after a polynomial fit (optional).
%   -Remove fast oscillations with a low-pass filter.
% 	-Normalize data as df/f or z-score according to settings.
% 	-Find peaks of spontaneous activity in low-pass filtered data.
%   -Figure 1:
%     -Raw signal, reference, and baseline model.
%     -Baseline corrected signal and reference.
%     -Motion correction.
%     -Normalization.
%     -Peak detection.
%   -Figure 2:
%     -Power spectrum for each epoch.
%   -Figure 3:
%     -Boxplot.
%   -Remaining figures:
%     -Peak triggers and event-triggers.
%     -Triggered average of peaks (if any) and event-triggers.
% 
% configuration is a struct with the following fields (defaults are used for missing fields):
%     conditionEpochs - Epochs for different conditions: {'epoch1', [start1, end1, start2, end2, ...], 'epoch2', ...}
%     artifactEpochs - Epochs to remove.
%     baselineEpochs - Epochs to include for baseline correction.
%     thresholdEpochs - Epochs to include for peak threshold calculation.
%     events - Event-triggered data; times at which a type of event occurs.
%     baselineLowpassFrequency - Frequency representative of baseline.
%     airPLS - Baseline correction for all data using airPLS (true, false, or airPLS inputs).
%     resamplingFrequency - Resampling frequency (Hz).
%     lowpassFrequency - Lowest frequency permitted in normalized signal.
%     peakSeparation - Minimum time separatio between two peaks.
%     peakWindow - Window to capture around each peak of spontaneous activity.
%     eventWindow - Window to capture around each event.
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
% 2022-06-30. Last modified.
classdef FPA < handle
    properties
        configuration
        
        warnings
        time
        frequency
        peakIds
        peakLabels
        eventIds
        eventLabels
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
        duration
        area
        peakThreshold
    end
    
    properties (Access = private)
        nConditions
        epochNames
        cleanId
        epochBounds
        epochLabels
        peakWindowTemplate
        eventWindowTemplate
        referenceProvided
        % Peaks detected including artifacts.
        peakIdsAll
        valleyIdsAll
        normalizedArea
        boxplotIds
        boxplotLabels
        eventNormalization
        
        % Settings for visualization.
        cmap = lines();
        zoomSettings = {0.99, 0.50};
        signalColor = [0.0000, 0.4470, 0.7410];
        referenceColor = [0.8500, 0.3250, 0.0980];
        alternateColor = [0, 0.6470, 0.9410];
        peaksLineColor = [0.4660, 0.6740, 0.1880];
        peaksMarkerColor = [1.0000, 0.0000, 0.0000];
        eventsMarkerColor = [0.0000, 0.0000, 0.0000];
        dashColor = [0.0000, 0.0000, 0.0000];
    end
    
    methods
        function obj = FPA(time, signal, reference, parameters)
            % Accumulate all warnings.
            obj.warnings = {};
            
            % The variable parameters is optional.
            if nargin < 4
                parameters = struct();
            end
            
            % Use defaults for missing parameters.
            defaults.conditionEpochs = {'Data', [-Inf, Inf]};
            defaults.artifactEpochs = [];
            defaults.baselineEpochs = [-Inf, Inf];
            defaults.events = [];
            defaults.resamplingFrequency = NaN;
            defaults.baselineLowpassFrequency = 0.1;
            defaults.airPLS = false;
            defaults.lowpassFrequency = 5.0;
            defaults.peakSeparation = 2.0;
            defaults.fitReference = true;
            defaults.f0 = @median;
            defaults.f1 = @mad;
            defaults.e0 = {@mean, [-1, 0]};
            defaults.e1 = 1;
            defaults.thresholdEpochs = NaN;
            defaults.threshold = {2.91, @mad, @median};
            defaults.peakWindow = 10.0;
            defaults.eventWindow = 10.0;
            
            % Override defaults with user parameters.
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
            
            if numel(configuration.peakWindow) == 1
                configuration.peakWindow = [-configuration.peakWindow, +configuration.peakWindow] / 2;
            end
            
            if numel(configuration.eventWindow) == 1
                configuration.eventWindow = [-configuration.eventWindow, +configuration.eventWindow] / 2;
            end

            % Resampling frequency defaults to the smallest between 100Hz and the source frequency.
            sourceFrequency = 1 / median(diff(time));
            if isnan(configuration.resamplingFrequency)
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
            elseif configuration.resamplingFrequency > sourceFrequency
                frequency = sourceFrequency;
                obj.warnings{end + 1} = warn('[resampling] Cannot resample to frequencies higher than the source frequency (%.2f Hz).', sourceFrequency);
            else
                frequency = sourceFrequency;
            end
            
            % Setup.
            nSamples = numel(time);
            nConditions = numel(configuration.conditionEpochs) / 2;
            epochNames = configuration.conditionEpochs(1:2:end);
            
            % Index of artifacts and non-artifacts.
            allIds = colon(1, nSamples)';
            artifactId = time2id(time, configuration.artifactEpochs);
            artifactFreeId = setdiff(allIds, artifactId);
            
            % Remove data from the edges as filtering will introduce artifacts here.
            % The exclusion window is the smallest between half a second and 5% of sample size.
            edgeWindow =  min(ceil(0.5 * frequency / configuration.lowpassFrequency), round(0.05 * nSamples));
            edgeId = union(artifactId, [1:edgeWindow, numel(time) - edgeWindow + 1:nSamples]');
            edgeFreeId = setdiff(allIds, edgeId);
            cleanId = intersect(artifactFreeId, edgeFreeId);
            
            % Define clean epochs for fitting and peak detection.
            baselineId = time2id(time, configuration.baselineEpochs);
            baselineId = intersect(baselineId, cleanId);
            
            % Threshold epochs defaults to everything.
            if isnan(configuration.thresholdEpochs)
                thresholdId = cleanId;
            else
                thresholdId = time2id(time, configuration.thresholdEpochs);
                thresholdId = intersect(cleanId, thresholdId);
            end
            
            % Replace artifacts with straight lines for modeling baseline.
            signalArtifactFree = signal;
            signalArtifactFree(artifactId) = interp1(artifactFreeId, signal(artifactFreeId), artifactId);
            referenceArtifactFree = reference;
            if referenceProvided
                referenceArtifactFree(artifactId) = interp1(artifactFreeId, reference(artifactFreeId), artifactId);
            end
            % Remove high-frequency oscillations to detect baseline (where indicated).
            signalArtifactFreeSmooth = signalArtifactFree;
            referenceArtifactFreeSmooth = referenceArtifactFree;
            if configuration.baselineLowpassFrequency <= frequency / 2
                baselineFilter = designfilt('lowpassiir', 'HalfPowerFrequency', configuration.baselineLowpassFrequency, 'SampleRate', frequency, 'DesignMethod', 'butter', 'FilterOrder', 12);
                signalArtifactFreeSmooth = filtfilt(baselineFilter, signalArtifactFree);
                if referenceProvided
                    referenceArtifactFreeSmooth = filtfilt(baselineFilter, referenceArtifactFree);
                end
            else
                obj.warnings{end + 1} = warn('[baseline-correction] Cannot lowpass to frequencies larger than half of the sampling frequency (%.2f Hz).', frequency / 2);
            end
            
            if useAirPLS
                % Model baseline with airPLS (everywhere).
                [~, signalBaseline] = airPLS(signalArtifactFreeSmooth', configuration.airPLS{:});
                signalBaseline = signalBaseline';
                signalCorrected = signalArtifactFree - signalBaseline;
                if referenceProvided
                    [~, referenceBaseline] = airPLS(referenceArtifactFreeSmooth', configuration.airPLS{:});
                    referenceBaseline = referenceBaseline';
                    referenceCorrected = referenceArtifactFree - referenceBaseline;
                else
                    referenceBaseline = zeros(size(signalCorrected));
                    referenceCorrected = zeros(size(signalCorrected));
                end
            else
                % Model baseline with an exponential decay at given epochs (where indicated).
                signalFit = fit(time(baselineId), signalArtifactFreeSmooth(baselineId), fittype('exp1'));
                signalBaseline = signalFit(time);
                signalCorrected = signalArtifactFree - signalBaseline;
                if referenceProvided
                    referenceFit = fit(time(baselineId), referenceArtifactFreeSmooth(baselineId), fittype('exp1'));
                    referenceBaseline = referenceFit(time);
                    referenceCorrected = referenceArtifactFree - referenceBaseline;
                else
                    referenceBaseline = zeros(size(signalCorrected));
                    referenceCorrected = zeros(size(signalCorrected));
                end
            end
            
            % Fit reference to signal (where indicated).
            if referenceProvided && configuration.fitReference
                r2sFit = fit(referenceCorrected(cleanId), signalCorrected(cleanId), fittype('poly1'), 'Robust', 'on');
                referenceCorrected = r2sFit.p1 * referenceCorrected + r2sFit.p2;
            end
            
            % Correct for movement artifacts.
            f = signalCorrected - referenceCorrected;
            
            % Low-pass filter.
            if configuration.lowpassFrequency == Inf
                fSmooth = f;
            else
                fFilter = designfilt('lowpassiir', 'HalfPowerFrequency', configuration.lowpassFrequency, 'SampleRate', frequency, 'DesignMethod', 'butter', 'FilterOrder', 12);
                fSmooth = filtfilt(fFilter, f);
            end
            
            % Normalize.
            [dff, f0, f1] = normalize(time, fSmooth, configuration.f0, configuration.f1);
            
            % Get peak threshold.
            state = warning('Query', 'signal:findpeaks:largeMinPeakHeight');
            warning('Off', 'signal:findpeaks:largeMinPeakHeight');
            peakThreshold = threshold(configuration.threshold, dff(thresholdId));
            
            if any(+dff >= peakThreshold)
                [~, peakIdsAll] = findpeaks(+dff, 'MinPeakHeight', peakThreshold, 'MinPeakDistance', configuration.peakSeparation * frequency);
            else
                peakIdsAll = [];
            end
            if any(-dff >= peakThreshold)
                [~, valleyIdsAll] = findpeaks(-dff, 'MinPeakHeight', peakThreshold, 'MinPeakDistance', configuration.peakSeparation * frequency);
            else
                valleyIdsAll = [];
            end
            warning(state.state, 'signal:findpeaks:largeMinPeakHeight');

            % Index template to apply around each peak and event.
            range = round(configuration.peakWindow * frequency);
            obj.peakWindowTemplate = range(1):range(2);
            range = round(configuration.eventWindow * frequency);
            obj.eventWindowTemplate = range(1):range(2);
            
            % Get indices for epochs.
            % Start and stop vector indices for all provided epochs.
            epochBounds = zeros(2, 0);
            % Numeric label corresponding to each epoch range.
            epochLabels = zeros(0, 1);
            % Vector index for each peak / event.
            peakIds = zeros(0, 1);
            eventIds = zeros(0, 1);
            % Numeric label corresponding to each peak / event.
            peakLabels = zeros(0, 1);
            eventLabels = zeros(0, 1);
            % Misc indexing / labeling.
            boxplotIds = zeros(0, 1);
            boxplotLabels = zeros(0, 1);
            peakCounts = zeros(nConditions, 1);
            valleyCounts = zeros(nConditions, 1);
            duration = zeros(nConditions, 1);
            area = zeros(nConditions, 1);
            
            % Get indices for time triggers.
            x = arrayfun(@(t) find(time >= t, 1, 'first'), configuration.events, 'UniformOutput', false);
            k = ~cellfun(@isempty, x);
            eventTimeIds = [x{k}];
            
            for c = 1:nConditions
                % Accumulate vector indices limited to conditions.
                [ids, bounds] = time2id(time, configuration.conditionEpochs{2 * c});
                
                % Start/stop-triggered data.
                n = numel(bounds) / 2;
                epochBounds = cat(2, epochBounds, bounds);
                epochLabels = cat(1, epochLabels, repmat(c, n, 1));
                
                % Event-triggered data.
                eventIdsEpoch = intersect(ids, eventTimeIds);
                n = numel(eventIdsEpoch);
                eventIds = cat(1, eventIds, eventIdsEpoch);
                eventLabels = cat(1, eventLabels, repmat(c, n, 1));
                
                % Peak-triggered data.
                peakIdsEpoch = intersect(ids, peakIdsAll);
                n = numel(peakIdsEpoch);
                peakIds = cat(1, peakIds, peakIdsEpoch);
                peakLabels = cat(1, peakLabels, repmat(c, n, 1));
                
                peakCounts(c) = sum(ismember(ids, peakIdsAll));
                valleyCounts(c) = sum(ismember(ids, valleyIdsAll));
                
                boxplotIds = cat(1, boxplotIds, ids);
                boxplotLabels = cat(1, boxplotLabels, repmat(c, numel(ids), 1));
                
                duration(c) = numel(ids) / frequency;
                area(c) = trapz(dff(ids)) / frequency;
            end

            % Normalize area according to epoch length.
            normalizedArea = area ./ duration;
            normalizedArea(duration == 0) = 0;
            
            obj.configuration = configuration;
            obj.time = time;
            obj.frequency = frequency;
            obj.eventNormalization = @(data) normalize(obj.eventWindowTemplate / obj.frequency, data, obj.configuration.e0, obj.configuration.e1);
            
            % Order depends on epoch definitions. Overlapping is possible and allowed.
            obj.peakIds = peakIds;
            obj.peakLabels = peakLabels;
            obj.eventIds = eventIds;
            obj.eventLabels = eventLabels;
            obj.epochBounds = epochBounds;
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
            obj.cleanId = edgeFreeId;
            obj.peakThreshold = peakThreshold;
            obj.peakIdsAll = peakIdsAll;
            obj.valleyIdsAll = valleyIdsAll;
            obj.boxplotIds = boxplotIds;
            obj.boxplotLabels = boxplotLabels;
            obj.normalizedArea = normalizedArea;
            obj.referenceProvided = referenceProvided;
        end
        
        function plot(obj)
            obj.plotTrace();
            %obj.plotPowerSpectrum();
            %obj.plotStatistics();
            obj.plotTrigger();
            obj.plotTriggerAverage();
            %obj.plotAUC();
        end
        
        function fig = plotTrace(obj)
            fig = figure('name', 'FPA: df/f');
            xlims = obj.time([1, end]);
            
            % Plot raw signal, reference, and baseline model.
            subplot(4, 1, 1);
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
            subplot(4, 1, 2);
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
            subplot(4, 1, 3);
            hold('all');
            yy = [obj.f(obj.cleanId); obj.fSmooth(obj.cleanId)];
            ylims = limits(yy, obj.zoomSettings{:});
            plotEpochs(obj.configuration.conditionEpochs, xlims, ylims, obj.cmap, false);
            plot(obj.time, obj.f, 'Color', obj.signalColor, 'DisplayName', 'f');
            plot(obj.time, obj.fSmooth, 'Color', obj.alternateColor, 'DisplayName', sprintf('f (<%.2fHz)', obj.configuration.lowpassFrequency));
            ylim(ylims);
            title('Motion correction');
            legend('show');
            
            % Plot normalization (e.g. df/f) and peak detection.
            subplot(4, 1, 4);
            hold('all');
            yy = obj.dff(obj.cleanId);
            ylims = limits(yy, obj.zoomSettings{:});
            epochs = obj.configuration.conditionEpochs;
            epochs(1:2:end) = arrayfun(@(e) sprintf('area:%.2f / %i peaks / %i valleys', obj.area(e), obj.peakCounts(e), obj.valleyCounts(e)), 1:obj.nConditions, 'UniformOutput', false);
            plotEpochs(epochs, xlims, ylims, obj.cmap, true);
            % Show event triggers if user provided such data, and the events exist within epochs.
            nEvents = numel(obj.eventIds);
            if nEvents > 0
                x = [obj.time(obj.eventIds), obj.time(obj.eventIds), NaN(nEvents, 1)]';
                y = repmat([ylims, NaN]', 1, nEvents);
                plot(x(:), y(:), 'Color', obj.eventsMarkerColor, 'LineStyle', '-', 'DisplayName', 'Events');
            end
            plot(obj.time, obj.dff, 'Color', obj.signalColor, 'DisplayName', 'df/f');
            plot(obj.time([1, end]), +obj.peakThreshold([1, 1]), 'Color', obj.dashColor, 'LineStyle', '--', 'DisplayName', sprintf('threshold:%.2f', obj.peakThreshold));
            plot(obj.time([1, end]), -obj.peakThreshold([1, 1]), 'Color', obj.dashColor, 'LineStyle', '--', 'HandleVisibility', 'off');
            plot(obj.time(obj.peakIdsAll), obj.dff(obj.peakIdsAll), 'Color', obj.peaksMarkerColor, 'LineStyle', 'none', 'Marker', 'o', 'DisplayName', 'Peaks');
            plot(obj.time(obj.valleyIdsAll), obj.dff(obj.valleyIdsAll), 'Color', obj.peaksMarkerColor, 'LineStyle', 'none', 'Marker', 'o', 'DisplayName', 'Peaks');
            ylim(ylims);
            title('Normalization and peak detection');
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
        
        function figs = plotTriggerAverage(obj)
            figs(1) = figure('name', 'FPA: Peak-triggered average');

            plotTriggerAverage(obj.dff, obj.peakIds, obj.peakLabels, obj.peakWindowTemplate, obj.frequency, @(data) data, obj.configuration.conditionEpochs(1:2:end), obj.cmap, 'No peaks');
            ylabel('df/f');
            title('Peak-triggered average');
            
            if numel(obj.eventIds) > 0
                figs(1) = figure('name', 'FPA: Event-triggered average');
                plotTriggerAverage(obj.dff, obj.eventIds, obj.eventLabels, obj.eventWindowTemplate, obj.frequency, obj.eventNormalization, obj.configuration.conditionEpochs(1:2:end), obj.cmap, 'No events');
                ylabel('df/f');
                title('Event-triggered average');
            end
        end
        
        function figs = plotTrigger(obj)
            name = 'Peak-trigger heatmap';
            figs(1) = figure('name', name);

            plotTriggerHeatmap(obj.dff, obj.peakIds, obj.peakLabels, obj.peakWindowTemplate, @(data) data, obj.frequency, obj.configuration.conditionEpochs(1:2:end), 'No peaks');
            annotation('textbox', [0, 0.95, 1, 0.05], 'string', name, 'LineStyle', 'none');
            
            if numel(obj.eventIds) > 0
                name = 'Event-trigger heatmap';
                figs(1) = figure('name', name);
                plotTriggerHeatmap(obj.dff, obj.eventIds, obj.eventLabels, obj.eventWindowTemplate, obj.eventNormalization, obj.frequency, obj.configuration.conditionEpochs(1:2:end), 'No events');
                annotation('textbox', [0, 0.95, 1, 0.05], 'string', name, 'LineStyle', 'none');
            end
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
            
            % Time vs f.
            % Rows represent increasing values of time with corresponding f values.
            output = fullfile(folder, sprintf('%s - f.csv', basename));
            fid = fopen(output, 'w');
            fprintf(fid, '# time, f\n');
            fprintf(fid, '%f, %f\n', [obj.time, obj.f]');
            fclose(fid);
            
            % Time vs dff.
            % Rows represent increasing values of time with corresponding dff values.
            output = fullfile(folder, sprintf('%s - dff.csv', basename));
            fid = fopen(output, 'w');
            fprintf(fid, '# time, dff\n');
            fprintf(fid, '%f, %f\n', [obj.time, obj.dff]');
            fclose(fid);
            
            % AUC.
            output = fullfile(folder, sprintf('%s - AUC.csv', basename));
            fid = fopen(output, 'w');
            fprintf(fid, '# conditionId, conditionName, area, duration, normalizedArea\n');
            data = [obj.epochNames(:), num2cell([colon(1, obj.nConditions)', obj.area, obj.duration, obj.normalizedArea])];
            data = data(:, [2, 1, 3, 4, 5])';
            fprintf(fid, '%i, %s, %f, %f, %f\n', data{:});
            fclose(fid);
            
            %  Compute fluorescence statistics for each epoch.
            output = fullfile(folder, sprintf('%s - epoch stats.csv', basename));
            fid = fopen(output, 'w');
            fprintf(fid, '# conditionId, conditionName, mean, median, max\n');
            startIds = obj.epochBounds(1:2:end);
            stopIds = obj.epochBounds(2:2:end);
            for c = 1:numel(startIds)
                x = obj.dff(startIds(c):stopIds(c));
                label = obj.epochLabels(c);
                fprintf(fid, '%i, %s, %f, %f, %f\n', obj.epochLabels(c), obj.epochNames{label}, mean(x), median(x), max(x));
            end
            fclose(fid);
            
            % All peak-triggered windows and their average with corresponding epoch label.
            output1 = fullfile(folder, sprintf('%s - peak-triggered.csv', basename));
            output2 = fullfile(folder, sprintf('%s - peak-triggered averaged.csv', basename));
            obj.saveEventTrigger(output1, output2, 'condition', obj.peakIds, obj.peakLabels, obj.peakWindowTemplate, onj.eventNormalization, obj.epochNames(obj.peakLabels));
            
            % All event-triggered windows and their average with corresponding epoch label.
            output1 = fullfile(folder, sprintf('%s - event-triggered.csv', basename));
            output2 = fullfile(folder, sprintf('%s - event-triggered averaged.csv', basename));
            obj.saveEventTrigger(output1, output2, 'condition', obj.eventIds, obj.eventLabels, obj.eventWindowTemplate, onj.eventNormalization, obj.epochNames(obj.eventLabels));
        end
    end
    
    methods (Access = private)
        function saveEventTrigger(obj, output1, output2, prefix, ids, labels, window, normalize, names)
            window = window(:);
            nTicks = numel(window);
            timeTicks = window / obj.frequency;
            timeHeader = strjoin(arrayfun(@(x) sprintf('%f', x), timeTicks, 'UniformOutput', false), ', ');
            
            % Filter out out-of-range traces.
            nSamples = numel(obj.dff);
            k = ids > min(window) & ids + max(window) < nSamples;
            labels = labels(k);
            ids = ids(k);
            ids = ids(:);
            labels = labels(:);
            if nargin == 7
                names = names(:);
                names = names(k);
                saveConditionName = true;
            else
                saveConditionName = false;
            end
            if numel(ids) >= 1
                windowIds = ids' + window;
                triggeredData = obj.dff;
                triggeredData = triggeredData(windowIds);
                triggeredData = reshape(triggeredData, nTicks, nTriggers);
                triggeredData = normalize(triggeredData);
                triggerTime = obj.time(ids);
                
                % All triggered windows with corresponding epoch label.
                % Order depends on epoch definitions. Overlapping is possible.
                % Rows represent a single trigger:
                % First column is the condition label of the trigger followed by the trace around each trigger, with each trigger at the center column (n / 2 + 1) labeled with "zero".
                fid = fopen(output1, 'w');
                if saveConditionName
                    fprintf(fid, ['# time, %1$sId, %1$sName, ', timeHeader, '\n'], prefix);
                    format = ['%f, %i, %s', repmat(', %f', 1, nTicks), '\n'];
                    data = [names, num2cell([triggerTime, labels, triggeredData])];
                    data = data(:, [2 3 1 4:end])';
                else
                    fprintf(fid, ['# time, %1$sId, ', timeHeader, '\n'], prefix);
                    format = ['%f, %i', repmat(', %f', 1, nTicks), '\n'];
                    data = num2cell([triggerTime, labels, triggeredData]);
                    data = data(:, :)';
                end
                fprintf(fid, format, data{:});
                fclose(fid);
                
                % Average of the above.
                uLabels = unique(labels);
                averages = zeros(0, size(triggeredData, 2));
                for u = 1:numel(uLabels)
                    label = uLabels(u);
                    uData = triggeredData(labels == label, :);
                    uData = reshape(uData, [], size(triggeredData, 2)); % !!
                    averages = cat(1, averages, mean(uData, 1));
                end
                
                fid = fopen(output2, 'w');
                if saveConditionName
                    fprintf(fid, ['# %1$sId, %1$sName, ', timeHeader, '\n'], prefix);
                    format = ['%i, %s', repmat(', %f', 1, nTicks), '\n'];
                    uTriggerNames = unique(names, 'stable');
                    data = [uTriggerNames, num2cell([uLabels, averages])];
                    data = data(:, [2, 1, 3:end])';
                    fprintf(fid, format, data{:});
                else
                    fprintf(fid, ['# %1$sId, ', timeHeader, '\n'], prefix);
                    format = ['%i', repmat(', %f', 1, nTicks), '\n'];
                    data = num2cell([uLabels, averages]);
                    data = data(:, :)';
                    fprintf(fid, format, data{:});
                end
                fclose(fid);
            end
        end
    end
end

function [data, f0, f1] = normalize(time, data, c0, c1)
    f0 = parseNormalization(time, data, c0);
    f1 = parseNormalization(time, data, c1);
    data = (data - f0) ./ f1;
end

function output = parseNormalization(time, data, parameters)
    % data is represented in a single column.
    % trials of data are represented in separate columns.
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
            nTicks = numel(time);
            window = parameters{2};
            window = min(round(window * frequency), nTicks);
            output = fcn(data, window, options{:});
        else
            % Produce a value from all data (or epochs).
            epochs = parameters{2};
            ids = time2id(time, epochs);
            output = fcn(data(ids, :));
        end
    elseif isa(parameters, 'function_handle')
        % Produce a value from all data (or epochs).
        fcn = parameters;
        epochs = [-Inf, Inf];
        ids = time2id(time, epochs);
        output = fcn(data(ids, :));
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
        epochs{2 * e}(epochs{2 * e} == -Inf) = xlims(1);
        epochs{2 * e}(epochs{2 * e} == +Inf) = xlims(2);
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

function plotTriggerHeatmap(data, ids, labels, window, normalize, frequency, names, message)
    % Filter out out-of-range traces.
    window = window(:);
    nSamples = numel(data);
    k = ids > min(window) & ids + max(window) < nSamples;
    labels = labels(k);
    ids = ids(k);
    nTicks = numel(window);
    timeTicks = window / frequency;
    
    nConditions = numel(names);
    [ni, nj] = squaredFactors(nConditions);
    
    if numel(ids) > 0
        allWindowIds = ids' + window;
        allTriggeredData = data(allWindowIds);
        clims = limits(allTriggeredData, 0.99, 0);
    end
    
    axs = cell(nConditions, 1);
    for c = 1:nConditions
        axs{c} = subplot(ni, nj, c);
        conditionIds = ids(labels == c);
        nTriggers = numel(conditionIds);
        if nTriggers > 0
            windowIds = conditionIds' + window;
            % Each column is an event.
            triggeredData = data(windowIds);
            % Make sure matrix is nr x nc, particularly when nr x 1.
            triggeredData = reshape(triggeredData, nTicks, nTriggers);
            triggeredData = normalize(triggeredData);
            imagesc('xData', timeTicks, 'yData', 1:nTriggers, 'cData', triggeredData', clims);
            yticks = get(gca(), 'YTick');
            yticks = yticks(round(yticks) == yticks);
            set(gca(), 'YTick', yticks);
        else
            text(0.5, 0.5, message, 'HorizontalAlignment', 'center');
        end
        title(names{c});
        axis('tight');
    end
    xlabel(axs{ni}, 'Time (s)');
    ylabel(axs{ni}, 'Trigger id');
    set(get(colorbar(), 'title'), 'string', 'df/f');
end

function plotTriggerAverage(data, ids, labels, window, frequency, normalize, names, colors, message)
    % Filter out out-of-range traces.
    window = window(:);
    nSamples = numel(data);
    k = ids > min(window) & ids + max(window) < nSamples;
    labels = labels(k);
    ids = ids(k);
    nTicks = numel(window);
    timeTicks = window / frequency;
    
    % Option 1: Difference of exponential functions.
    %model = @(p, x) p(1) + p(2) * (exp(-x / p(3)) - exp(-x / p(4)));
    %mask = true(size(time));
    %initial = @(y, duration) [median(y), std(y), duration / 2, duration / 2];
    
    % Option 2: Exponential function.
    model = @(p, x) p(1) + p(2) * exp(-x / p(3));
    mask = timeTicks >= 0;
    initial = @(y, duration) [median(y), std(y), duration / 2];

    duration = max(timeTicks) - min(timeTicks);
    options = optimoptions('lsqcurvefit', 'Algorithm', 'levenberg-marquardt', 'Display', 'off');
    
    if numel(ids) > 0
        hold('all');
        nConditions = numel(names);
        for c = 1:nConditions
            conditionIds = ids(labels == c);
            nTriggers = numel(conditionIds);
            if nTriggers > 0
                windowIds = conditionIds' + window;
                % Each column is an event.
                triggeredData = data(windowIds);
                % Make sure matrix is nr x nc, particularly when nr x 1.
                triggeredData = reshape(triggeredData, nTicks, nTriggers);
                triggeredData = normalize(triggeredData);
                av = mean(triggeredData, 2);
                sd = std(triggeredData, [], 2);
                
                x = timeTicks(mask) - min(timeTicks(mask));
                y = av(mask);
                parameters = lsqcurvefit(model, initial(y, duration), x(:), y(:), [], [], options);
                tauFall = parameters(end);
                
                % Plot.
                plot(timeTicks, av, 'Color', colors(c, :), 'HandleVisibility', 'off');
                if parameters(2) <= 0
                    tauFallText = sprintf('<Inf>');
                else
                    tauFallText = sprintf('%.2fs', tauFall);
                    plot(timeTicks(mask), model(parameters, timeTicks(mask)), 'Color', colors(c, :), 'LineStyle', '--', 'HandleVisibility', 'off');
                end
                sem = sd / sqrt(nTriggers);
                sem0 = sem(ceil(nTriggers / 2));
                vertices = [timeTicks, av + sem / 2];
                vertices = cat(1, vertices, [flipud(timeTicks), flipud(av - sem / 2)]);
                faces = 1:2 * numel(window);
                label = sprintf('%s (SEM=%.2f, tau=%s, n = %i)', names{c}, sem0, tauFallText, nTriggers);
                patch('Faces', faces, 'Vertices', vertices, 'FaceColor', colors(c, :), 'EdgeColor', 'none', 'FaceAlpha', 0.10, 'DisplayName', label);
            end
        end
    else
        text(0.5, 0.5, message, 'HorizontalAlignment', 'center');
    end
    plot(NaN(2, 1), NaN(2, 1), 'Color', 'k', 'LineStyle', '--', 'DisplayName', 'Exponential fit');
    legend('show');
    xlabel('Time (s)');
    axis('tight');
end