% FPA - Fiber Photometry Analysis.
% 
% fpa = FPA(time, signal, reference, configuration);
% Subtract photobleaching, correct for motion artifacts, 
% normalize, filter, detect peaks of spontaneous activity in user defined
% epochs, display trigger averages of user defined events.
% 
% Units for time and frequency are seconds and hertz respectively.
% 
% See examples and documentation online for detailed analysis steps and
% default parameters: https://github.com/leomol/FPA
% 
% FPA methods:
% FPA - Class constructor
% 
% FPA methods to plot results:
% plotAUC - Plot area under the curve
% plotStatistics - Make a boxplot
% plotTrace - Plot pre-processing steps
% plotPowerSpectrum - Plot power spectrum
% plotEvents - Plot event-triggered average
% plotEventHeatmap - Plot event-triggered heat map
% plotPeaks - Plot peak-triggered average
% plotPeakHeatmap - Plot peak-triggered heat map
% 
% FPA methods to export preprocessed data:
% exportF - Export pre-processed f and normalized f
% 
% FPA methods to export triggered data:
% exportEvents - Export event-triggered traces
% exportEventAverage - Export event-triggered average
% exportPeaks - Export peak-triggered traces
% exportPeakAverage - Export peak-triggered average
% 
% FPA methods to export epoch delimited data:
% exportStatistics - Export AUC, sum, mean, median, min and max for normalized f and epoch duration and peak counts.
% 
% FPA methods to customize the configuration:
% adaptive - Get a baseline using the AirPLS function using default parameters
% findPeaks - Find peaks in data
% fit - Fit curve to trace according to the fitting type
% fitReferenceToSignal - Fit reference to signal using non-negative, least-squares fit, at the given epochs
% get - Get trace according to data choice and epoch range
% ids - Get indices relative to timeResampled for the given epochs
% interpolate - Interpolate trace linearly within epochs
% lowpass - Low-pass filter trace according to frequency
% normalize - Normalize data according to parameters f0 and f1
% resample - Resample signal and reference according to frequency
% 
% FPA static methods:
% Adaptive - Get a baseline using the AirPLS function using default parameters
% Defaults - Get default configuration
% Fit - Fit curve to trace according to the fitting type
% Ids - Get indices relative to time for the given epochs
% Interpolate - Interpolate trace linearly within epochs
% Lowpass - Low-pass filter trace according to frequency
% Normalize - Normalize data according to parameters f0 and f1
% 
% 2019-02-01. Leonardo Molina.
% 2023-12-07. Last modified.
classdef FPA < handle
    properties (Access = public)
        % time - Raw time
        time
        % signal - Raw signal
        signal
        % reference - Raw reference
        reference
        % timeResampled - Time after resampling
        timeResampled
        % signalResampled - Signal after resampling
        signalResampled
        % referenceResampled - Reference after resampling
        referenceResampled
        % signalTrimmed - Signal after removing artifacts
        signalTrimmed
        % referenceArtifactFree - Reference after removing artifacts
        referenceTrimmed
        % signalSmoothed - Smoothed signal to model photobleaching
        signalSmoothed
        % referenceSmoothed - Smoothed reference to model photobleaching
        referenceSmoothed
        % signalModeled - Baseline modeled from signal
        signalModeled
        % referenceModeled - Baseline modeled from reference
        referenceModeled
        % signalCorrected - Signal after subtracting baseline
        signalCorrected
        % referenceCorrected - Reference after subtracting baseline
        referenceCorrected
        % signalStandardized - Standardized signal
        signalStandardized
        % signalStandardized - Standardized reference
        referenceStandardized
        % referenceFitted - Reference fitted to signal
        referenceFitted
        % f - F
        f
        % fSmoothed - Smoothed f
        fSmoothed
        % fNormalized - df/f
        fNormalized
        % f0 - f0
        f0
        % f1 - f1
        f1
    end

    properties (SetAccess = private)
        % epochs - Cell array containing epoch name and epoch ranges.
        epochs
        
        % resampleData - Function that returns resampled time, signal, and reference
        resampleData
        % trimSignal - Function that removes artifacts from the signal
        trimSignal
        % trimReference - Function that removes artifacts from the reference
        trimReference
        % smoothSignal - Function that smooths fpa.signalTrimmed, for the purposes of modeling the baseline
        smoothSignal
        % smoothReference - Function that smooths fpa.referenceTrimmed, for the purposes of modeling the baseline
        smoothReference
        % modelSignal - Function that models the baseline by fitting a curve to fpa.signalSmoothed
        modelSignal
        % modelReference - Function that models the baseline by fitting a curve to fpa.referenceSmoothed
        modelReference
        % correctSignal - Function that removes the baseline from the signal
        correctSignal
        % correctReference - Function that removes the baseline from the reference
        correctReference
        % standardizeSignal - Function that standardizes the signal
        standardizeSignal
        % standardizeReference - Function that standardizes the reference
        standardizeReference
        % fitReference - Function that rescales the reference so that it is comparable to the signal
        fitReference
        % getF - Function that returns the motion- and bleaching corrected data
        getF
        % smoothF - unction that smooths f
        smoothF
        % normalizeF - Function that normalizes f
        normalizeF
        % normalizeEvents - Function that normalizes each event trace
        normalizeEvents
        % normalizePeaks - Function that normalizes each peak trace
        normalizePeaks
        % getPeaks - Function to detect peaks for peak-triggered averages
        getPeaks
        
        % duration - Duration of each epoch
        duration
        % area - Area under the curve for each epoch
        area
        % normalizedArea - Area under the curve normalized by the epoch duration
        normalizedArea
        
        % peakIds - Index of each peak detected
        peakIds
        % peakLabels - Labels given to peaks found within the same epoch
        peakLabels
        % peakCounts - Number of peaks within each epoch
        peakCounts
        
        % frequency - Resampling frequency
        frequency
    end
    
    properties (Access = private)
        % nEpochs - Number of epoch definitions or conditions
        nEpochs
        % epochNames - Name of each epoch definition extracted from epochs
        epochNames
        % epochRanges - Epoch time ranges extracted from epochs
        epochRanges
        % epochStartIds - Start indices for each period and epoch of epochRanges
        epochStartIds
        % epochStopIds - Stop indices for each period and epoch of epochRanges
        epochStopIds
        % epochStartLabels - Labels given to periods belonging to the same epochRange
        epochStartLabels
        % referenceProvided - Whether a reference was provided
        referenceProvided
        % epochIds - Indices corresponding to all epoch definitions
        epochIds
        % epochLabels - Labels given to indices found within the same epoch
        epochLabels
        
        % Settings for visualization.
        cmap = lines();
        zoomSettings = {99.75, 0.50};
        colors = struct();
    end
    
    methods
        function obj = FPA(time, signal, reference, parameters)
            % FPA(time, signal, reference, parameters)
            % Preprocess signal and reference according to options defined in the parameters.
            
            % The variable parameters is optional.
            if nargin < 4
                parameters = struct();
            end
            
            obj.colors.signal = [0.0000, 0.4470, 0.7410];
            obj.colors.reference = [0.8500, 0.3250, 0.0980];
            obj.colors.referenceFitted = [0.0000, 0.0000, 0.0000];
            obj.colors.fSmoothed = [0, 0.6470, 0.9410];
            obj.colors.peaks = [1.0000, 0.0000, 0.0000];
            obj.colors.threshold = [0.0000, 0.0000, 0.0000];
            
            % Get defaults.
            defaults = FPA.Defaults();
            
            % Override defaults with user parameters.
            expectedParameters = fieldnames(defaults);
            providedParameters = fieldnames(parameters);
            remainingParameters = setdiff(expectedParameters, providedParameters);
            for i = 1:numel(remainingParameters)
                name = remainingParameters{i};
                obj.(name) = defaults.(name);
            end
            terminate = false;
            for i = 1:numel(providedParameters)
                name = providedParameters{i};
                if ismember(name, expectedParameters)
                    obj.(name) = parameters.(name);
                else
                    terminate = true;
                    warn('[parsing] "%s" is not a valid parameter.', name);
                end
            end
            if terminate
                error('See previous warnings.');
            end
            
            % Setup.
            obj.epochNames = obj.epochs(1:2:end);
            obj.epochRanges = obj.epochs(2:2:end);
            obj.nEpochs = numel(obj.epochNames);
            
            % Turn data into columns.
            obj.time = time(:);
            obj.signal = signal(:);
            obj.reference = reference(:);
            
            % Reference may be empty.
            obj.referenceProvided = ~isempty(reference);
            
            % Resample data to target frequency.
            if isEnabled(obj.resampleData)
                [obj.timeResampled, obj.signalResampled, obj.referenceResampled] = obj.resampleData(obj);
            else
                [obj.timeResampled, obj.signalResampled, obj.referenceResampled] = deal(obj.time, obj.signal, obj.reference);
            end
            obj.frequency = getFrequency(obj.timeResampled);
            
            % Replace artifacts with straight lines for modeling baseline.
            if isEnabled(obj.trimSignal)
                obj.signalTrimmed = obj.trimSignal(obj, obj.timeResampled, obj.signalResampled);
            else
                obj.signalTrimmed = obj.signalResampled;
            end
            if obj.referenceProvided
                if isEnabled(obj.trimReference)
                    obj.referenceTrimmed = obj.trimReference(obj, obj.timeResampled, obj.referenceResampled);
                else
                    obj.referenceTrimmed = obj.referenceResampled;
                end
            end
            
            % Remove high-frequency oscillations to detect baseline (where indicated).
            if isEnabled(obj.smoothSignal)
                obj.signalSmoothed = obj.smoothSignal(obj, obj.timeResampled, obj.signalTrimmed);
            else
                obj.signalSmoothed = obj.signalTrimmed;
            end
            if obj.referenceProvided
                if isEnabled(obj.smoothReference)
                    obj.referenceSmoothed = obj.smoothReference(obj, obj.timeResampled, obj.referenceTrimmed);
                else
                    obj.referenceSmoothed = obj.referenceTrimmed;
                end
            end
            
            % Detect baseline according to model.
            if isEnabled(obj.modelSignal)
                obj.signalModeled = obj.modelSignal(obj, obj.timeResampled, obj.signalSmoothed);
            else
                obj.signalModeled = zeros(size(obj.signalSmoothed));
            end
            if obj.referenceProvided
                if isEnabled(obj.modelReference)
                    obj.referenceModeled = obj.modelReference(obj, obj.timeResampled, obj.referenceSmoothed);
                else
                    obj.referenceModeled = zeros(size(obj.referenceSmoothed));
                end
            end
            
            % Subtract baseline.
            if isEnabled(obj.correctSignal)
                obj.signalCorrected = obj.correctSignal(obj);
            else
                obj.signalCorrected = obj.signalTrimmed;
            end
            if obj.referenceProvided
                if isEnabled(obj.correctReference)
                    obj.referenceCorrected = obj.correctReference(obj);
                else
                    obj.referenceCorrected = obj.referenceTrimmed;
                end
            end

            % Standardize signal and reference.
            if isEnabled(obj.standardizeSignal)
                obj.signalStandardized = obj.standardizeSignal(obj);
            else
                obj.signalStandardized = obj.signalCorrected;
            end
            if obj.referenceProvided
                if isEnabled(obj.standardizeReference)
                    obj.referenceStandardized = obj.standardizeReference(obj);
                else
                    obj.referenceStandardized = obj.referenceCorrected;
                end
            end
            
            % Fit reference to signal.
            if obj.referenceProvided
                if isEnabled(obj.fitReference)
                    obj.referenceFitted = obj.fitReference(obj);
                else
                    obj.referenceFitted = obj.referenceStandardized;
                end
            end
            
            % Unfiltered f.
            if isEnabled(obj.getF)
                obj.f = obj.getF(obj);
            else
                obj.f = obj.signalStandardized;
            end
            
            % Filtered, motion corrected.
            if isEnabled(obj.smoothF)
                obj.fSmoothed = obj.smoothF(obj, obj.timeResampled, obj.f);
            else
                obj.fSmoothed = obj.f;
            end
            
            % Normalize.
            if isEnabled(obj.normalizeF)
                [obj.fNormalized, obj.f0, obj.f1] = obj.normalizeF(obj, obj.timeResampled, obj.fSmoothed);
            else
                [obj.fNormalized, obj.f0, obj.f1] = deal(obj.f, NaN, NaN);
            end
            
            if isEnabled(obj.getPeaks)
                obj.peakIds = obj.getPeaks(obj);
            else
                obj.peakIds = [];
            end
            
            % Start and stop vector indices for all provided epochs.
            obj.epochStartIds = zeros(0, 1);
            obj.epochStopIds = zeros(0, 1);
            % Numeric label corresponding to each epoch range.
            obj.epochStartLabels = zeros(0, 1);
            % Numeric label corresponding to each peak.
            obj.peakLabels = zeros(size(obj.peakIds));
            % Misc indexing / labeling.
            obj.epochIds = zeros(0, 1);
            obj.epochLabels = zeros(0, 1);
            obj.peakCounts = zeros(obj.nEpochs, 1);
            obj.duration = zeros(obj.nEpochs, 1);
            obj.area = zeros(obj.nEpochs, 1);
            
            for c = 1:obj.nEpochs
                % Accumulate vector indices limited to conditions.
                [k, bounds] = obj.ids(obj.epochRanges{c});
                startIds = bounds(1, :);
                stopIds = bounds(2, :);
                % Start/stop-triggered data.
                n = numel(bounds) / 2;
                obj.epochStartIds = cat(1, obj.epochStartIds, startIds(:));
                obj.epochStopIds = cat(1, obj.epochStopIds, stopIds(:));
                obj.epochStartLabels = cat(1, obj.epochStartLabels, repmat(c, n, 1));
                
                % Peak-triggered data.
                obj.peakLabels(ismember(obj.peakIds, k)) = c;
                obj.peakCounts(c) = sum(ismember(k, obj.peakIds));
                
                obj.epochIds = cat(1, obj.epochIds, k);
                obj.epochLabels = cat(1, obj.epochLabels, repmat(c, numel(k), 1));
                
                obj.duration(c) = numel(k) / obj.frequency;
                obj.area(c) = trapz(obj.fNormalized(k)) / obj.frequency;
            end
            
            % Normalize area according to epoch length.
            obj.normalizedArea = obj.area ./ obj.duration;
            obj.normalizedArea(obj.duration == 0) = 0;
        end

        function data = adaptive(~, data, varargin)
            % FPA.adaptive(data, lambda, order, wep, p, itermax)
            % Model a baseline in data using the airPLS algorithm with default parameters.
            data = FPA.Adaptive(data, varargin{:});
        end

        function peakIds = findPeaks(obj, type, separation, peakThreshold)
            % peakIds = FPA.findPeaks(type, separation, peakThreshold)
            % Shortcut to MATLAB's findpeaks.
            
            state = warning('Query', 'signal:findpeaks:largeMinPeakHeight');
            warning('Off', 'signal:findpeaks:largeMinPeakHeight');
            options = {'MinPeakHeight', 'MinPeakProminence'};
            [~, peakIds] = findpeaks(obj.fNormalized, options{strcmpi({'height', 'prominence'}, type)}, peakThreshold, 'MinPeakDistance', separation * obj.frequency);
            warning(state.state, 'signal:findpeaks:largeMinPeakHeight');
        end

        function data = fit(~, time, data, fitType, epochs)
            % data = FPA.fit(time, data, fitType, epochs)
            % Fit a curve of fitType to a given trace, at the given epochs.
            
            data = FPA.Fit(time, data, fitType, epochs);
        end
        
        function data1 = fitReferenceToSignal(obj, targetFrequency, epochs)
            % data = FPA.fitReferenceToSignal(targetFrequency, epochs)
            % Fit reference to the smoothed signal using a non-negative, least-squares fit, at the given epochs.

            if nargin < 2
                targetFrequency = 0.1;
            end
            if nargin < 3
                epochs = [-Inf, Inf];
            end
            
            % Fit on smoothed data.
            k = obj.ids(epochs);
            data1 = FPA.Lowpass(obj.timeResampled, obj.referenceStandardized, targetFrequency);
            data2 = FPA.Lowpass(obj.timeResampled, obj.signalStandardized, targetFrequency);
            data1 = obj.referenceStandardized * lsqnonneg(data1(k), data2(k));
        end

        function data = get(obj, traceName, epochs)
            % data = FPA.get(traceName, epochs)
            % Get a trace at the given epoch ranges.
            % 
            % Example:
            %   Get reference data at the given epochs:
            %   epochs = [0, 10, 100, 200];
            %   Option 1:
            %     fpa.reference(fpa.ids(epochs))
            %   Option 2:
            %     fpa.get('reference', epoch)
            
            data = obj.(traceName);
            k = obj.ids(epochs);
            data = data(k);
        end
        
        function [k, limits] = ids(obj, epochs)
            % [ids, bounds] = FPA.ids(epochs)
            % Returns the index in time where time is enclosed by the given epochs
            
            [k, limits] = FPA.Ids(obj.timeResampled, epochs);
        end

        function data = interpolate(~, time, data, epochs)
            % data = FPA.interpolate(time, data, epochs)
            % Remove artifacts from the given trace at the given epochs.
            
            data = FPA.Interpolate(time, data, epochs);
        end
        
        function data = lowpass(~, time, data, targetFrequency)
            % data = FPA.lowpass(time, data, frequency)
            % Low-pass filter the given trace at the given frequency.
            
            data = FPA.Lowpass(time, data, targetFrequency);
        end
        
        function [data, f0, f1] = normalize(~, time, data, f0, f1)
            % [data, f0, f1] = FPA.normalize(time, data, f0, f1)
            % Normalize data as (f - f0) / f1 where f0 and f1 are one of
            %   1) A scalar (0, 1, ...)
            %   2) An array of the same shape as time ([0, 0, 0, 0, ...], ...)
            %   3) A function handle returning an scalar (@mean, @median, @std, @mad)
            %   4) A cell array indicating a handle to a non-moving function and an epoch array ({@mean, [0, 100, 300, 350]}, ...).
            %   5) A cell array indicating a handle to a moving function and a window size ({@movmean, 60}, ...).
            
            [data, f0, f1] = FPA.Normalize(time, data, f0, f1);
        end

        function [time, signal, reference] = resample(obj, frequency)
            % [time, signal, reference] = FPA.resample(frequency)
            % Resample data to target frequency.
            
            sourceFrequency = getFrequency(obj.time);
            if frequency < sourceFrequency
                targetFrequency = frequency;
                % Express frequency as a ratio p/q.
                [p, q] = rat(targetFrequency / sourceFrequency);
                % Resample: interpolate every p/q/f, upsample by p, filter, downsample by q.
                [signal, time2] = resample(obj.signal, obj.time, targetFrequency, p, q);
                k = time2 <= obj.time(end);
                signal = signal(k);
                if ~isempty(obj.reference)
                    reference = resample(obj.reference, obj.time, targetFrequency, p, q);
                    reference = reference(k);
                end
                time = time2(k);
            else
                time = obj.time;
                signal = obj.signal;
                reference = obj.reference;
            end
        end
        
        function fig = plotAUC(obj, normalized)
            % fig = FPA.plotAUC(normalized)
            % Plot area under the curve.
            % Area is normalized when normalized is set to true (default).

            normalized = nargin < 2 || normalized;
            
            % Plot normalized area under the curve.
            name = 'Normalized area under the curve';
            fig = figure('name', ['FPA' name]);
            ax.auc = axes();
            if normalized
                data = obj.normalizedArea;
            else
                data = obj.area;
            end
            bar(1:obj.nEpochs, data);
            set(ax.auc, 'XTickLabel', obj.epochNames);
            xtickangle(45);
            if normalized
                ylabel('Norm. area');
                title('dff/f - normalized AUC');
            else
                ylabel('Area');
                title('dff/f - AUC');
            end
        end

        function fig = plotPowerSpectrum(obj, overlap, window)
            % fig = FPA.plotPowerSpectrum(overlap, window)
            % Plot power spectrum of df/f.
            % The calculation is done in traces of a size given by window (default is 10s).
            % Set overlap to true to plot all epochs on a single plot.
            
            if nargin < 2
                overlap = true;
            end
            if nargin < 3
                window = 10;
            end

            name = 'Power spectrum';
            fig = figure('name', ['FPA: ' name]);
            axs = cell(1, obj.nEpochs);
            segmentLength = round(window * obj.frequency);
            for c = 1:obj.nEpochs
                if overlap
                    axs{c} = gca;
                else
                    axs{c} = subplot(obj.nEpochs, 1, c);
                end
                hold('all');
                epochName = obj.epochNames{c};
                epochRange = obj.epochRanges{c};
                k = obj.ids(epochRange);
                n = numel(k);
                if n > 2
                    [power, frequencies, ci] = pwelch(obj.fNormalized(k), ones(segmentLength, 1), 0, segmentLength, obj.frequency, 'power', 'ConfidenceLevel', 0.95);
                    vertices = [frequencies, 10 * log10(ci(:, 2))];
                    vertices = cat(1, vertices, flipud([frequencies, 10 * log10(ci(:, 1))]));
                    faces = 1:2 * numel(frequencies);
                    patch('Faces', faces, 'Vertices', vertices, 'FaceColor', obj.cmap(c, :), 'EdgeColor', 'none', 'FaceAlpha', 0.10, 'HandleVisibility', 'off');
                    plot(frequencies, 10 * log10(power), 'Color', obj.cmap(c, :), 'DisplayName', epochName);
                end
                if overlap
                    title('Power spectrum');
                else
                    subName = sprintf('%s - Power spectrum', epochName);
                    title(subName);
                end
                axis('tight');
            end
            if overlap
                legend('show');
            end
            axs = [axs{:}];
            ylims = get(axs, 'ylim');
            ylims = cat(1, ylims{:});
            ylims = [min(ylims(:, 1)), max(ylims(:, 2))];
            set(axs, 'ylim', ylims);
            ylabel('Power (dB)');
            xlabel('Frequency (Hz)');
            linkaxes(axs, 'x');
        end
        
        function fig = plotStatistics(obj)
            % fig = FPA.plotStatistics()
            % Box plot of traces for each epoch.

            % Boxplot of fNormalized.
            name = 'Boxplot';
            fig = figure('name', ['FPA: ' name]);
            % Not all epochs may be available.
            boxplotNames = obj.epochNames(unique(obj.epochLabels));
            boxplot(obj.fNormalized(obj.epochIds), obj.epochLabels, 'Labels', boxplotNames);
            hold('all');
            ylims = ylim();
            for c = 1:obj.nEpochs
                k = obj.ids(obj.epochRanges{c});
                n = numel(k);
                if n > 2
                    epochStatLabel = sprintf('\nmean:%.2f\nstd:%.2f', obj.normalizedArea(c), std(obj.fNormalized(k)));
                    text(c, ylims(2), epochStatLabel, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'Top');
                end
            end
            plot(obj.normalizedArea, 'k.', 'DisplayName', 'Mean');
            ylabel('df/f');
            xtickangle(45);
            title('Stats on df/f traces for each epoch');
        end
        
        function fig = plotTrace(obj)
            % fig = FPA.plotTrace()
            % Plot traces for each stage of the analysis.

            name = 'df/f';
            fig = figure('name', ['FPA:' name]);
            xlims = obj.timeResampled([1, end]);
            
            % Plot raw signal, reference, and photobleaching model.
            subplot(5, 1, 1);
            hold('all');
            yy = [obj.signalResampled(:); obj.referenceResampled(:); obj.signalModeled(:)];
            ylims = limits(yy, obj.zoomSettings{:});
            plotEpochs(obj.epochNames, obj.epochRanges, xlims, ylims, obj.cmap, true);
            if obj.referenceProvided
                plot(obj.timeResampled, obj.referenceResampled, 'Color', obj.colors.reference, 'DisplayName', 'Resampled reference');
            end
            plot(obj.timeResampled, obj.signalResampled, 'Color', obj.colors.signal, 'DisplayName', 'Resampled signal');
            if obj.referenceProvided
                plot(obj.timeResampled, obj.referenceModeled, 'Color', obj.colors.threshold, 'LineStyle', '--', 'DisplayName', 'Photobleaching', 'HandleVisibility', 'off');
            end
            plot(obj.timeResampled, obj.signalModeled, 'Color', obj.colors.threshold, 'LineStyle', '--', 'DisplayName', 'Photobleaching');
            ylim(ylims);
            title('Raw data');
            legend('show');
            
            % Plot bleaching corrected signal and reference.
            subplot(5, 1, 2);
            hold('all');
            yy = [obj.signalCorrected(:); obj.referenceCorrected(:)];
            ylims = limits(yy, obj.zoomSettings{:});
            plotEpochs(obj.epochNames, obj.epochRanges, xlims, ylims, obj.cmap, false);
            if obj.referenceProvided
                plot(obj.timeResampled, obj.referenceCorrected, 'Color', obj.colors.reference, 'DisplayName', 'Corrected reference');
            end
            plot(obj.timeResampled, obj.signalCorrected, 'Color', obj.colors.signal, 'DisplayName', 'Corrected signal');
            ylim(ylims);
            title('Photobleaching correction');
            legend('show');

            % Plot standardize signal and reference.
            subplot(5, 1, 3);
            hold('all');
            yy = [obj.signalStandardized(:); obj.referenceStandardized(:)];
            ylims = limits(yy, obj.zoomSettings{:});
            plotEpochs(obj.epochNames, obj.epochRanges, xlims, ylims, obj.cmap, false);
            if obj.referenceProvided
                plot(obj.timeResampled, obj.referenceStandardized, 'Color', obj.colors.reference, 'DisplayName', 'Standardized reference');
            end
            plot(obj.timeResampled, obj.signalStandardized, 'Color', obj.colors.signal, 'DisplayName', 'Standardized signal');
            if obj.referenceProvided
                if isEnabled(obj.fitReference)
                    plot(obj.timeResampled, obj.referenceFitted, 'Color', obj.colors.referenceFitted, 'DisplayName', 'Fitted reference');
                end
            end
            ylim(ylims);
            title('Standardization');
            legend('show');
            
            % Plot motion correction (f and lowpass f).
            subplot(5, 1, 4);
            hold('all');
            yy = [obj.f; obj.fSmoothed];
            ylims = limits(yy, obj.zoomSettings{:});
            plotEpochs(obj.epochNames, obj.epochRanges, xlims, ylims, obj.cmap, false);
            plot(obj.timeResampled, obj.f, 'Color', obj.colors.threshold, 'DisplayName', 'f');
            plot(obj.timeResampled, obj.fSmoothed, 'Color', obj.colors.signal, 'DisplayName', 'f smooth');
            ylim(ylims);
            title('Motion correction');
            legend('show');
            
            % Plot normalization (e.g. df/f) and peak detection.
            subplot(5, 1, 5);
            hold('all');
            yy = obj.fNormalized;
            ylims = limits(yy, obj.zoomSettings{:});
            tmpNames = arrayfun(@(e) sprintf('area:%.2f / %i peaks', obj.area(e), obj.peakCounts(e)), 1:obj.nEpochs, 'UniformOutput', false);
            plotEpochs(tmpNames, obj.epochRanges, xlims, ylims, obj.cmap, true);
            plot(obj.timeResampled, obj.fNormalized, 'Color', obj.colors.signal, 'DisplayName', 'df/f');
            ids = obj.peakIds(obj.peakLabels > 0);
            plot(obj.timeResampled(ids), obj.fNormalized(ids), 'Color', obj.colors.peaks, 'LineStyle', 'none', 'Marker', 'o', 'DisplayName', 'Peaks');
            ylim(ylims);
            title('Normalization and peak detection');
            legend('show');
            
            % Move axes together.
            axs = findobj(fig, 'type', 'axes');
            linkaxes(axs, 'x');
            xlim(axs(1), [obj.timeResampled(1), obj.timeResampled(end)]);
            set(axs(2:end), 'XTick', []);
            xlabel('Time (s)');
            ylabel('df/f');
        end
        
        function fig = plotEvents(obj, eventTimes, window)
            % fig = FPA.plotEvents(eventTimes, window)
            % Plot event-triggered average around a window (default is 5s).

            if nargin < 3
                window = 5;
            end
            name = 'Event-triggered average';
            fig = figure('name', ['FPA: ' name]);
            [eventIds, eventLabels] = obj.getEventData(eventTimes);
            if numel(eventIds) > 0
                template = obj.parseWindow(window);
                plotTriggerAverage(obj.fNormalized, eventIds, eventLabels, template, obj.frequency, @(time, data) obj.normalizeEvents(obj, time, data), obj.epochNames, obj.cmap, 'No events'); % !!
                ylabel('df/f');
                title(name);
            end
        end
        
        function fig = plotEventHeatmap(obj, eventTimes, window)
            % fig = FPA.plotEventHeatmap(eventTimes, window)
            % Plot event-triggered heatmap around a window (default is 5s).

            if nargin < 3
                window = 5;
            end
            name = 'Event-trigger heatmap';
            fig = figure('name', ['FPA: ' name]);
            [eventIds, eventLabels] = obj.getEventData(eventTimes);
            if numel(eventIds) > 0
                template = obj.parseWindow(window);
                plotTriggerHeatmap(obj.fNormalized, eventIds, eventLabels, template, @(time, data) obj.normalizeEvents(obj, time, data), obj.frequency, obj.epochNames, 'No events'); % !!
                annotation('textbox', [0, 0.95, 1, 0.05], 'string', name, 'LineStyle', 'none');
            end
        end
        
        function fig = plotPeaks(obj, window)
            % fig = FPA.plotPeaks(window)
            % Plot peak-triggered average around a window (default is 5s).

            if nargin < 2
                window = 5;
            end
            name = 'Peak-triggered average';
            fig = figure('name', ['FPA: ' name]);
            template = obj.parseWindow(window);
            k = obj.peakLabels > 0;
            plotTriggerAverage(obj.fNormalized, obj.peakIds(k), obj.peakLabels(k), template, obj.frequency, @(time, data) obj.normalizePeaks(obj, time, data), obj.epochNames, obj.cmap, 'No peaks'); % !! obj.epochNames(obj.peakLabels(k)) ?
            ylabel('df/f');
            title(name);
        end
        
        function fig = plotPeakHeatmap(obj, window)
            % fig = FPA.plotPeakHeatmap(window)
            % Plot peak-triggered heatmap around a window (default is 5s).

            if nargin < 2
                window = 5;
            end
            name = 'Peak-trigger heatmap';
            fig = figure('name', ['FPA: ' name]);
            template = obj.parseWindow(window);
            k = obj.peakLabels > 0;
            plotTriggerHeatmap(obj.fNormalized, obj.peakIds(k), obj.peakLabels(k), template, @(time, data) obj.normalizePeaks(obj, time, data), obj.frequency, obj.epochNames, 'No peaks'); % !!
            annotation('textbox', [0, 0.95, 1, 0.05], 'string', name, 'LineStyle', 'none');
        end
        
        function exportF(obj, filename)
            % FPA.exportF(filename)
            % Export time vs f  (both normalized and unnormalized).
            % Rows represent increasing values of time with corresponding dff values.
            
            fid = fopen(filename, 'w');
            fprintf(fid, 'time, f, dff\n');
            fprintf(fid, '%f, %f, %f\n', [obj.timeResampled, obj.f, obj.fNormalized]');
            fclose(fid);
        end
        
        function exportStatistics(obj, filename)
            % FPA.exportStatistics(filename)
            % Export AUC, sum, mean, median, min and max for normalized f and epoch duration and peak counts.

            fid = fopen(filename, 'w');
            fprintf(fid, 'conditionId, conditionName, duration, area, sum, mean, median, min, max, peakCount\n');
            for c = 1:numel(obj.epochStartIds)
                x = obj.fNormalized(obj.epochStartIds(c):obj.epochStopIds(c));
                label = obj.epochStartLabels(c);
                fprintf(fid, '%i, %s, %f, %f, %f, %f, %f, %f, %f, %i\n', obj.epochStartLabels(c), obj.epochNames{label}, obj.duration(c), obj.area(c), sum(x), mean(x), median(x), min(x), max(x), obj.peakCounts(c));
            end
            fclose(fid);
        end
        
        function exportEvents(obj, filename, eventTimes, window, asColumns)
            % FPA.exportEvents(filename, eventTimes, window, asColumns)
            % Export event-triggered traces with corresponding epoch name and label.
            
            if nargin < 4
                window = 5;
            end
            
            if nargin < 5
                asColumns = true;
            end
            template = obj.parseWindow(window);
            [eventIds, eventLabels] = obj.getEventData(eventTimes);
            obj.saveTriggers(filename, asColumns, false, eventIds, eventLabels, template, @(time, data) obj.normalizeEvents(obj, time, data), obj.epochNames(eventLabels));
        end
        
        function exportEventAverage(obj, filename, eventTimes, window, asColumns)
            % FPA.exportEventAverage(filename, eventTimes, window, asColumns)
            % Export event-triggered average with corresponding epoch name and label.

            if nargin < 4
                window = 5;
            end
            if nargin < 5
                asColumns = true;
            end

            template = obj.parseWindow(window);
            [eventIds, eventLabels] = obj.getEventData(eventTimes);
            obj.saveTriggers(filename, asColumns, true, eventIds, eventLabels, template, @(time, data) obj.normalizeEvents(obj, time, data), obj.epochNames(eventLabels));
        end
        
        function exportPeaks(obj, filename, window, asColumns)
            % FPA.exportPeaks(filename, window, asColumns)
            % All peak-triggered traces with corresponding epoch name and label.

            if nargin < 3
                window = 5;
            end
            if nargin < 4
                asColumns = true;
            end
            
            template = obj.parseWindow(window);
            k = obj.peakLabels > 0;
            obj.saveTriggers(filename, asColumns, false, obj.peakIds(k), obj.peakLabels(k), template, @(time, data) obj.normalizePeaks(obj, time, data), obj.epochNames(obj.peakLabels(k)));
        end
        
        function exportPeakAverage(obj, filename, window, asColumns)
            % FPA.exportPeakAverage(filename, window, asColumns)
            % All peak-triggered average with corresponding epoch name and label.

            if nargin < 3
                window = 5;
            end
            if nargin < 4
                asColumns = true;
            end
            
            template = obj.parseWindow(window);
            k = obj.peakLabels > 0;
            obj.saveTriggers(filename, asColumns, true, obj.peakIds(k), obj.peakLabels(k), template, @(time, data) obj.normalizePeaks(obj, time, data), obj.epochNames(obj.peakLabels(k)));
        end
    end
    
    methods (Access = private)
        function [eventIds, eventLabels] = getEventData(obj, eventTimes)
            % [eventIds, eventLabels] = FPA.getEventData(eventTimes)
            % Get eventIds and eventLabels given eventTimes.

            eventIds = zeros(0, 1);
            eventLabels = zeros(0, 1);
            
            % Get indices for time triggers.
            x = arrayfun(@(t) find(obj.timeResampled >= t, 1, 'first'), eventTimes, 'UniformOutput', false);
            r = ~cellfun(@isempty, x);
            eventTimeIds = [x{r}];
            
            for c = 1:obj.nEpochs
                % Accumulate vector indices limited to conditions.
                k = obj.ids(obj.epochRanges{c});
                
                % Event-triggered data.
                eventIdsEpoch = intersect(k, eventTimeIds);
                n = numel(eventIdsEpoch);
                eventIds = cat(1, eventIds, eventIdsEpoch);
                eventLabels = cat(1, eventLabels, repmat(c, n, 1));
            end
        end
        
        function template = parseWindow(obj, window)
            % template = FPA.parseWindow(window)
            % Get an indexing template to apply around each peak and event.

            n = numel(window);
            if n == 1
                template = -round(0.5 * window * obj.frequency):round(0.5 * window * obj.frequency);
            elseif n == 2
                template = round(window(1) * obj.frequency):round(window(2) * obj.frequency);
            else
                template = round(window * obj.frequency);
            end
        end
        
        function saveTriggers(obj, output, asColumns, saveAverage, ids, labels, window, normalization, names)
            % FPA.saveTriggers(output, asColumns, saveAverage, ids, labels, window, normalization, names)
            % Save peak or event triggers to disk.

            window = window(:);
            nTicks = numel(window);
            timeTicks = window / obj.frequency;
            timeHeader = strjoin(arrayfun(@(x) sprintf('%f', x), timeTicks, 'UniformOutput', false), ', ');
            
            % Filter out out-of-range traces.
            nSamples = numel(obj.fNormalized);
            k = ids + min(window) > 0 & ids + max(window) < nSamples;
            labels = labels(k);
            ids = ids(k);
            ids = ids(:);
            labels = labels(:);
            names = names(:);
            names = names(k);
            uniqueNames = unique(names, 'stable');
            nUniqueNames = numel(uniqueNames);
            nTriggers = numel(ids);
            if nTriggers >= 1
                windowIds = ids' + window;
                triggeredData = obj.fNormalized;
                triggeredData = triggeredData(windowIds);
                triggeredData = reshape(triggeredData, nTicks, nTriggers);
                triggeredData = normalization(timeTicks, triggeredData);
                
                % Average of the above.
                uLabels = unique(labels);
                n = numel(uLabels);
                averages = zeros(nTicks, n);
                for u = 1:n
                    label = uLabels(u);
                    uData = triggeredData(:, labels == label);
                    % Make sure it has nTicks rows even when nTriggers = 1.
                    uData = reshape(uData, nTicks, []);
                    uData = mean(uData, 2);
                    averages(:, u) = uData;
                end
                
                % All triggered windows with corresponding epoch label.
                % Order depends on epoch definitions. Overlapping is possible.
                % Rows represent a single trigger:
                % First column is the condition label of the trigger followed by the trace around each trigger, with each trigger at the center column (n / 2 + 1) labeled with "zero".
                
                fid = fopen(output, 'w');
                if saveAverage
                    if asColumns
                        % time, epochName1, epochName2, ...
                        % -0.3,              x,              x, ...
                        % -0.2,              x,              x, ...
                        % -0.1,              x,              x, ...
                        fprintf(fid, ['time, ', strjoin(uniqueNames, ', '), '\n']);
                        format = ['%f', repmat(', %f', 1, nUniqueNames), '\n'];
                        data = num2cell([timeTicks, averages])';
                    else
                        % epochId, epochName, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3
                        %       1,         A,    x,    x,    x,   x,   x,   x,   x
                        %       2,         B,    x,    x,    x,   x,   x,   x,   x
                        %       3,         C,    x,    x,    x,   x,   x,   x,   x
                        fprintf(fid, ['epochId, epochName, ', timeHeader, '\n']);
                        format = ['%i, %s', repmat(', %f', 1, nTicks), '\n'];
                        data = [num2cell(uLabels)'; uniqueNames'; num2cell(averages)];
                    end
                else
                    if asColumns
                        %     , epochName1, ... epochName2, epochName2, ... 
                        % time,       0.00,           5.91,      9.91, ... 
                        % -0.3,          x,                         x, ...
                        % -0.2,          x,                         x, ...
                        % -0.1,          x,                         x, ...
                        fprintf(fid, ['    , ', strjoin(names, ', '), '\n']);
                        triggerTimeHeader = strjoin(arrayfun(@(x) sprintf('%f', x), obj.timeResampled(ids), 'UniformOutput', false), ', ');
                        fprintf(fid, ['time, ', triggerTimeHeader, '\n']);
                        format = ['%f', repmat(', %f', 1, nTriggers), '\n'];
                        data = num2cell([timeTicks, triggeredData])';
                    else
                        % time, epochId, epochName, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3
                        % 0.00,       1,         A,    x,    x,    x,   x,   x,   x,   x
                        % 5.91,       2,         B,    x,    x,    x,   x,   x,   x,   x
                        % 9.91,       2,         B,    x,    x,    x,   x,   x,   x,   x
                        fprintf(fid, ['time, epochId, epochName, ', timeHeader, '\n']);
                        format = ['%f, %i, %s', repmat(', %f', 1, nTicks), '\n'];
                        data = [num2cell(obj.timeResampled(ids))'; num2cell(labels)'; names'; num2cell(triggeredData)];
                    end
                end
                fprintf(fid, format, data{:});
                fclose(fid);
            end
        end
    end
    
    properties (Constant)
        % version - FPA version
        version = '2.0.4'

        % defaults - Configuration defaults.
        defaults = FPA.Defaults();
    end

    methods (Static)
        function config = Defaults()
            % FPA.Defaults()
            % Default parameters for preprocessing.
            
            config = struct();
            config.epochs = {'Data', [-Inf, Inf]};
            config.resampleData = @(fpa) fpa.resample(100);
            config.trimSignal = @(fpa, time, data) fpa.interpolate(time, data, []);
            config.trimReference = @(fpa, time, data) fpa.interpolate(time, data, []);
            config.smoothSignal = @(fpa, time, data) fpa.lowpass(time, data, 0.1);
            config.smoothReference = @(fpa, time, data) fpa.lowpass(time, data, 0.1);
            config.modelSignal = @(fpa, time, data) fpa.fit(time, data, 'exp1', [-Inf, Inf]);
            config.modelReference = @(fpa, time, data) fpa.fit(time, data, 'exp1', [-Inf, Inf]);
            config.correctSignal = @(fpa) fpa.signalTrimmed - fpa.signalModeled;
            config.correctReference = @(fpa) fpa.referenceTrimmed - fpa.referenceModeled;
            config.standardizeSignal = @(fpa) zscore(fpa.signalCorrected);
            config.standardizeReference = @(fpa) zscore(fpa.referenceCorrected);
            config.fitReference = @(fpa) fpa.fitReferenceToSignal(0.1, [-Inf, Inf]);
            config.getF = @(fpa) fpa.signalStandardized - fpa.referenceFitted;
            config.smoothF = @(fpa, time, data) fpa.lowpass(time, data, 10);
            config.normalizeF = @(fpa, time, data) fpa.normalize(time, data, @median, @mad);
            config.normalizeEvents = @(fpa, time, data) fpa.normalize(time, data, {@mean, [-Inf, 0]}, 1);
            config.normalizePeaks = @(fpa, time, data) fpa.normalize(time, data, {@mean, [-Inf, 0]}, 1);
            config.getPeaks = @(fpa) fpa.findPeaks('prominence', 0.5, median(fpa.fNormalized) + 2.91 * mad(fpa.fNormalized));
        end
        
        function data = Adaptive(data, varargin)
            % FPA.Adaptive(data, lambda, order, wep, p, itermax)
            % Model a baseline in data using the airPLS algorithm with default parameters.
            
            parameters = {5e9, 2, 0.1, 0.5, 50};
            [parameters{1:numel(varargin)}] = deal(varargin{:});
            if ~isempty(data)
                [~, data] = FPA.airPLS(data(:)', parameters{:});
                data = data(:);
            end
        end

        function data = Fit(time, data, fitType, epochs)
            % data = FPA.Fit(time, data, fitType, epochs)
            % Fit a curve of fitType to a given trace, at the given epochs.
            
            k = FPA.Ids(time, epochs);
            fitModel = fit(time(k), data(k), fittype(fitType));
            data = fitModel(time);
        end
        
        function [ids, limits] = Ids(time, epochs)
            % [ids, bounds] = FPA.Ids(time, epochs)
            % Returns the index in time where time is enclosed by the given epochs

            % (from1, to1, from2, to2, ...).
            % 2019-02-01. Leonardo Molina.
            % 2022-06-06. Last modified.

            epochs = epochs(:);
            % timeLimits = zeros(2, nEpochs);
            a = arrayfun(@(t) find(time >= t, 1, 'first'), epochs(1:2:end), 'UniformOutput', false);
            b = arrayfun(@(t) find(time <= t, 1, 'last'), epochs(2:2:end), 'UniformOutput', false);
            k = cellfun(@isempty, a) | cellfun(@isempty, b);
            a(k) = [];
            b(k) = [];
            nEpochs = numel(a);
            limits = zeros(2, nEpochs);
            limits(1:2:end) = [a{:}];
            limits(2:2:end) = [b{:}];
            % When "last" can't find anything ahead of "first", force a single time point.
            k = diff(limits, [], 1) < 0;
            limits(2, k) = limits(1, k);
            ids = arrayfun(@(e) colon(limits(1, e), limits(2, e))', 1:nEpochs, 'UniformOutput', false);
            ids = cat(1, ids{:});
        end

        function data = Interpolate(time, data, epochs)
            % data = FPA.Interpolate(time, data, epochs)
            % Remove artifacts from the given trace at the given epochs.

            if ~isempty(epochs)
                % Index of artifacts and non-artifacts.
                artifactId = FPA.ids(time, epochs);
                artifactFreeId = setdiff(colon(1, numel(time))', artifactId);
                data(artifactId) = interp1(artifactFreeId, data(artifactFreeId), artifactId);
            end
        end
        
        function data = Lowpass(time, data, targetFrequency)
            % data = FPA.Lowpass(time, data, frequency)
            % Low-pass filter the given trace at the given frequency.
            
            order = 12;
            sourceFrequency = getFrequency(time);
            filter = designfilt('lowpassiir', 'HalfPowerFrequency', targetFrequency, 'SampleRate', sourceFrequency, 'DesignMethod', 'butter', 'FilterOrder', order);
            data = filtfilt(filter, data);
        end
        
        function [data, f0, f1] = Normalize(time, data, parametersF0, parametersF1)
            % [data, f0, f1] = FPA.Normalize(time, data, f0, f1)
            % Normalize data as (data - f0) / f1
            f0 = parseNormalizationParameter(time, data, parametersF0);
            f1 = parseNormalizationParameter(time, data, parametersF1);
            data = (data - f0) ./ f1;
        end
    end
end

function output = parseNormalizationParameter(time, data, parameters)
    % output = parseNormalizationParameter(time, data, parameters)
    % parameters: value, [value1, value2, ...], @mean, {@mean, [start1, end1, start2, end2, ...]}, {@movemean, windowSize}, @median, ...
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
            frequency = getFrequency(time);
            nTicks = numel(time);
            window = parameters{2};
            window = min(round(window * frequency), nTicks);
            output = fcn(data, window, options{:});
        else
            % Produce a value from all data (or epochs).
            epochRange = parameters{2};
            ids = FPA.Ids(time, epochRange);
            output = fcn(data(ids, :));
        end
    elseif isa(parameters, 'function_handle')
        % Produce a value from all data (or epochs).
        fcn = parameters;
        epochRange = [-Inf, Inf];
        ids = FPA.Ids(time, epochRange);
        output = fcn(data(ids, :));
    else
        output = parameters;
    end
end

function plotEpochs(epochNames, epochRanges, xlims, ylims, cmap, show)
    for e = 1:numel(epochRanges)
        epochName = epochNames{e};
        epochRange = epochRanges{e};
        epochRange(epochRange == -Inf) = xlims(1);
        epochRange(epochRange == +Inf) = xlims(2); % !!
        [faces, vertices] = patchEpochs(epochRange, ylims(1), ylims(2));
        vertices(vertices == -Inf) = xlims(1);
        vertices(vertices == +Inf) = xlims(2); % !!
        if show
            patch('Faces', faces, 'Vertices', vertices, 'FaceColor', cmap(e, :), 'EdgeColor', 'none', 'FaceAlpha', 0.50, 'DisplayName', sprintf('%s', epochName)); % !!
        else
            patch('Faces', faces, 'Vertices', vertices, 'FaceColor', cmap(e, :), 'EdgeColor', 'none', 'FaceAlpha', 0.50, 'HandleVisibility', 'off');
        end
    end
end

function warn(format, varargin)
    fprintf(2, '[%s] %s\n', mfilename(), sprintf(format, varargin{:}));
end

function ylims = limits(x, percentile, grow)
    x = x(:);
    ylims = [prctile(x, 100 - percentile), prctile(x, percentile)];
    delta = diff(ylims) * grow;
    ylims = [ylims(1) - delta, ylims(2) + delta];
    ylims = [max(min(x), ylims(1)), min(max(x), ylims(2))];
end

function plotTriggerHeatmap(data, ids, labels, window, normalization, frequency, names, message)
    % Filter out out-of-range traces.
    window = window(:);
    nSamples = numel(data);
    k = ids + min(window) > 0 & ids + max(window) < nSamples;
    labels = labels(k);
    ids = ids(k);
    nTicks = numel(window);
    timeTicks = window / frequency;
    
    nEpochs = numel(names);
    [ni, nj] = squaredFactors(nEpochs);
    
    if numel(ids) > 0
        allWindowIds = ids' + window;
        allTriggeredData = data(allWindowIds);
        clims = limits(allTriggeredData, 99.75, 0);
    end
    
    axs = cell(nEpochs, 1);
    for c = 1:nEpochs
        axs{c} = subplot(ni, nj, c);
        conditionIds = ids(labels == c);
        nTriggers = numel(conditionIds);
        if nTriggers > 0
            windowIds = conditionIds' + window;
            % Each column is an event.
            triggeredData = data(windowIds);
            % Make sure matrix is nr x nc, particularly when nr x 1.
            triggeredData = reshape(triggeredData, nTicks, nTriggers);
            triggeredData = normalization(timeTicks, triggeredData);
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

function plotTriggerAverage(data, ids, labels, window, frequency, normalization, names, colors, message)
    % Filter out out-of-range traces.
    window = window(:);
    nSamples = numel(data);
    r = ids + min(window) > 0 & ids + max(window) < nSamples;
    labels = labels(r);
    ids = ids(r);
    nTicks = numel(window);
    timeTicks = window / frequency;
    if numel(ids) > 0
        hold('all');
        nEpochs = numel(names);
        for c = 1:nEpochs
            conditionIds = ids(labels == c);
            nTriggers = numel(conditionIds);
            if nTriggers > 0
                windowIds = conditionIds' + window;
                % Each column is an event.
                triggeredData = data(windowIds);
                % Make sure matrix is nr x nc, particularly when nr x 1.
                triggeredData = reshape(triggeredData, nTicks, nTriggers);
                triggeredData = normalization(timeTicks, triggeredData);
                av = mean(triggeredData, 2);
                sd = std(triggeredData, [], 2);
                % Plot.
                plot(timeTicks, av, 'Color', colors(c, :), 'HandleVisibility', 'off');
                sem = sd / sqrt(nTriggers);
                sem0 = sem(ceil(numel(sem) / 2));
                vertices = [timeTicks, av + sem / 2];
                vertices = cat(1, vertices, flipud([timeTicks, av - sem / 2]));
                faces = 1:2 * nTicks;
                label = sprintf('%s (SEM=%.2f, n=%i)', names{c}, sem0, nTriggers);
                patch('Faces', faces, 'Vertices', vertices, 'FaceColor', colors(c, :), 'EdgeColor', 'none', 'FaceAlpha', 0.10, 'DisplayName', label);
            end
        end
    else
        text(0.5, 0.5, message, 'HorizontalAlignment', 'center');
    end
    legend('show');
    xlabel('Time (s)');
    axis('tight');
end

function result = isEnabled(value)
    result = isa(value, 'function_handle');
end

% [faces, vertices] = patchEpochs(epochs, minimum, maximum)
% Create vectors faces and vertices from vector epochs (from1, to1, from2, to2, ...)
% to feed to MATLAB's function patch. The y-limits of the "drawing" are
% given by minimum and maximum.
% 2019-02-01. Leonardo Molina.
% 2019-08-30. Last modified.
function [faces, vertices] = patchEpochs(epochs, minimum, maximum)
    epochs = epochs(:);
    nEpochs = numel(epochs);
    vertices = zeros(2, 2 * nEpochs);
    for e = 1:2:nEpochs
        range = 4 * (e - 1) + (1:8);
        e1 = e;
        e2 = e + 1;
        vertices(range) = [epochs(e1), minimum, epochs(e1), maximum, epochs(e2), maximum, epochs(e2), minimum];
    end
    vertices = vertices';
    faces = reshape(1:2 * nEpochs, 4, nEpochs / 2)';
end

% [ni, nj] = squaredFactors(number, landscape)
% Return two closest divisors of a number.
% 2021-03-15. Leonardo Molina.
% 2021-03-15. Last modified.
function [ni, nj] = squaredFactors(number, landscape)
    landscape = nargin == 1 || landscape;
    number = round(number);
    x = 1:number;
    rr = x(~(rem(number, x)));
    cc = number ./ rr;
    [~, k] = mink(abs(rr - cc), 1);
    if landscape == (cc(k) < rr(k))
        ni = cc(k);
        nj = rr(k);
    else
        ni = rr(k);
        nj = cc(k);
    end
end

function frequency = getFrequency(time)
    frequency = 1 / median(diff(time));
end