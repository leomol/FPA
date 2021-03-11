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
% 2021-03-11. Last modified.
classdef FPA < handle
    properties
        configuration
        
        warnings
        figures
        time
        frequency
        peakIds
        peakLabels
        epochIds
        epochLabels
        
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
            defaults.plot = true;
            
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

            % Settings for visualization.
            percentile = 0.99;
            grow = 0.50;

            % Sampling frequency defaults.
            sourceFrequency = 1 / median(diff(time));
            if ~ismember('resamplingFrequency', parametersNames)
                configuration.resamplingFrequency = min(100, sourceFrequency);
            end

            % Plot: true | false | cell array of choices.
            if islogical(configuration.plot)
                if configuration.plot
                    configuration.plot = {'trace', 'power', 'stats', 'trigger', 'AUC'};
                else
                    configuration.plot = {};
                end
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
            valleyThreshold = -peakThreshold;

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
            windowTemplate = -halfWindow:halfWindow;

            % Index epochs.
            epochIds = zeros(2, 0);
            epochLabels = zeros(0, 1);
            peakIds = zeros(0, 1);
            peakLabels = zeros(0, 1);
            boxplotIds = zeros(0, 1);
            boxplotLabels = zeros(0, 1);
            peakCount = zeros(nConditions, 1);
            valleyCount = zeros(nConditions, 1);
            area = zeros(nConditions, 1);
            duration = zeros(nConditions, 1);
            normalizedArea = zeros(nConditions, 1);

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

                peakCount(c) = sum(ismember(ids, uniquePeakIds));
                valleyCount(c) = sum(ismember(ids, uniqueValleyIds));

                duration(c) = numel(ids) / frequency;
                area(c) = trapz(dff(ids)) / frequency;
                if numel(ids) > 0
                    normalizedArea(c) = area(c) / duration(c);
                end
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
            obj.figures = [];

            anyMatch = @(choices, pattern) any(~cellfun(@isempty, regexp(choices, pattern, 'start', 'once')));
            if anyMatch(configuration.plot, '\<trace\>')
                obj.figures(end + 1) = figure('name', 'FPA: df/f');

                % Plot raw signal, reference, and baseline model.
                ax.raw = subplot(5, 1, 1);
                ax.raw.XTick = [];
                hold(ax.raw, 'all');
                yy = [signal(:); reference(:); signalBaseline(:)];
                ylims = limits(yy, percentile, grow);
                plotEpochs(configuration.conditionEpochs, xlims, ylims, cmap, true);
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
                plotEpochs(configuration.conditionEpochs, xlims, ylims, cmap, false);
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
                plotEpochs(configuration.conditionEpochs, xlims, ylims, cmap, false);
                plot(ax.f, time, f, 'Color', signalColor, 'DisplayName', 'f');
                plot(ax.f, time, fSmooth, 'Color', alternateColor, 'DisplayName', sprintf('f (<%.2fHz)', configuration.lowpassFrequency));
                ylim(ax.f, ylims);
                title(ax.f, 'Motion correction');
                legend(ax.f, 'show');

                % Plot normalization (e.g. df/f).
                ax.filtered = subplot(5, 1, 4);
                ax.filtered.XTick = [];
                hold(ax.filtered, 'all');
                yy = [dff(cleanIds); peaksSmooth(cleanIds)];
                ylims = limits(yy, percentile, grow);
                epochs = configuration.conditionEpochs;
                epochs(1:2:end) = arrayfun(@(e) sprintf('area:%.2f', area(e)), 1:nConditions, 'UniformOutput', false);
                plotEpochs(epochs, xlims, ylims, cmap, true);
                plot(ax.filtered, time, dff, 'Color', signalColor, 'DisplayName', 'df/f');
                plot(ax.filtered, time, peaksSmooth, 'Color', peaksLineColor, 'DisplayName', sprintf('df/f (<%.2fHz)', configuration.peaksLowpassFrequency));
                ylim(ax.filtered, ylims);
                title(ax.filtered, 'Normalization');
                legend(ax.filtered, 'show');

                % Plot peak detection.
                ax.processed = subplot(5, 1, 5);
                hold(ax.processed, 'all');
                yy = peaksSmooth(cleanIds);
                ylims = limits(yy, percentile, grow);
                epochs = configuration.conditionEpochs;
                epochs(1:2:end) = arrayfun(@(e) sprintf('%i peaks / %i valleys', peakCount(e), valleyCount(e)), 1:nConditions, 'UniformOutput', false);
                plotEpochs(epochs, xlims, ylims, cmap, true);
                plot(ax.processed, time, peaksSmooth, 'Color', peaksLineColor, 'DisplayName', sprintf('df/f (<%.2fHz)', configuration.peaksLowpassFrequency));
                plot(ax.processed, time([1, end]), peakThreshold([1, 1]), 'Color', dashColor, 'LineStyle', '--', 'DisplayName', sprintf('threshold:%.2f', peakThreshold));
                plot(ax.processed, time([1, end]), valleyThreshold([1, 1]), 'Color', dashColor, 'LineStyle', '--', 'HandleVisibility', 'off');
                plot(ax.processed, time(uniquePeakIds), peaksSmooth(uniquePeakIds), 'Color', peaksMarkerColor, 'LineStyle', 'none', 'Marker', 'o', 'HandleVisibility', 'off');
                plot(ax.processed, time(uniqueValleyIds), peaksSmooth(uniqueValleyIds), 'Color', peaksMarkerColor, 'LineStyle', 'none', 'Marker', 'o', 'HandleVisibility', 'off');
                ylim(ylims);
                title(ax.processed, 'Peak detection');
                legend(ax.processed, 'show');

                % Move axes together.
                linkaxes(findobj(gcf(), 'type', 'axes'), 'x');
                xlim(ax.raw, [time(1), time(end)]);

                xlabel('Time (s)');
                ylabel('df/f');
            end

            if anyMatch(configuration.plot, '\<power\>')
                % Plot power spectrum.
                obj.figures(end + 1) = figure('name', 'FPA: Power spectrum');
                axs = cell(1, nConditions);
                for c = 1:nConditions
                    axs{c} = subplot(nConditions, 1, c);
                    epochName = configuration.conditionEpochs{2 * c - 1};
                    ids = time2id(time, configuration.conditionEpochs{2 * c});
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

            if anyMatch(configuration.plot, '\<stats\>')
                % Boxplot of dff.
                obj.figures(end + 1) = figure('name', 'FPA: Boxplot');
                ax.boxplot = axes();
                % Not all epochs may be available.
                boxplotNames = epochNames(unique(boxplotLabels));
                boxplot(dff(boxplotIds), boxplotLabels, 'Labels', boxplotNames);
                hold('all');
                ylims = ylim();
                for c = 1:nConditions
                    ids = time2id(time, configuration.conditionEpochs{2 * c});
                    n = numel(ids);
                    if n > 2
                        epochStatLabel = sprintf('\nmean:%.2f\nstd:%.2f', normalizedArea(c), std(dff(ids)));
                        text(c, ylims(2), epochStatLabel, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'Top');
                    end
                end
                plot(normalizedArea, 'k.', 'DisplayName', 'Mean');
                ylabel('df/f');
                xtickangle(45);
                title('Stats on df/f traces for each condition');
            end

            if anyMatch(configuration.plot, '\<trigger\>')
                obj.figures(end + 1) = figure('name', 'FPA: Peak-triggered average');
                plotTriggerAverage(dff, peakIds, peakLabels, windowTemplate, frequency, cmap, configuration.conditionEpochs(1:2:end));
                ylabel('df/f');
                title('Peak-triggered average');

                obj.figures(end + 1) = figure('name', 'FPA: Epoch start-triggered average');
                plotTriggerAverage(dff, epochIds(1:2:end)', epochLabels, windowTemplate, frequency, cmap, configuration.conditionEpochs(1:2:end));
                ylabel('df/f');
                title('Epoch start-triggered average');

                obj.figures(end + 1) = figure('name', 'FPA: Epoch stop-triggered average');
                plotTriggerAverage(dff, epochIds(2:2:end)', epochLabels, windowTemplate, frequency, cmap, configuration.conditionEpochs(1:2:end));
                ylabel('df/f');
                title('Epoch stop-triggered average');
            end

            if anyMatch(configuration.plot, '\<AUC\>')
                % Plot normalized area under the curve.
                nConditions = numel(configuration.conditionEpochs) / 2;
                obj.figures(end + 1) = figure('name', 'FPA: Normalized area under the curve');
                ax.auc = axes();
                bar(1:nConditions, normalizedArea);
                set(ax.auc, 'XTickLabel', epochNames);
                xtickangle(45);
                title('dff/f - normalized AUC');
            end
            
            obj.configuration = configuration;
            obj.time = time;
            obj.frequency = frequency;
            
            % Order depends on epoch definitions. Overlapping is possible and allowed.
            obj.peakIds = peakIds;
            obj.peakLabels = peakLabels;
            obj.epochIds = epochIds;
            obj.epochLabels = epochLabels;

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
        end
        
        function export(obj, prefix)
            prefix = regexprep(prefix, '/$', '');
            [folder, basename] = fileparts(prefix);
            
            % Save data for post-processing.
            nConditions = numel(obj.configuration.conditionEpochs) / 2;

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
            fprintf(fid, '# condition, area, duration\n');
            fprintf(fid, '%i, %.4f, %d\n', [(1:nConditions)', obj.area, obj.duration]');
            fclose(fid);

            % All peak-triggered windows and their average with corresponding epoch label.
            output1 = fullfile(folder, sprintf('%s - peak-triggered.csv', basename));
            output2 = fullfile(folder, sprintf('%s - peak-triggered averaged.csv', basename));
            obj.saveEventTrigger(output1, output2, obj.peakIds, obj.peakLabels);

            % Same as above for start-triggered windows.
            output1 = fullfile(folder, sprintf('%s - start-triggered.csv', basename));
            output2 = fullfile(folder, sprintf('%s - start-triggered averaged.csv', basename));
            obj.saveEventTrigger(output1, output2, obj.epochIds(1:2:end)', obj.epochLabels);

            % Same as above for stop-triggered windows.
            output1 = fullfile(folder, sprintf('%s - stop-triggered.csv', basename));
            output2 = fullfile(folder, sprintf('%s - stop-triggered averaged.csv', basename));
            obj.saveEventTrigger(output1, output2, obj.epochIds(2:2:end)', obj.epochLabels);
        end
        
        function saveEventTrigger(obj, output1, output2, triggerIds, triggerLabels)
            [triggeredWindow, halfWindow] = forceOdd(obj.configuration.triggeredWindow * obj.frequency);
            window = -halfWindow:halfWindow;
            traceTime = window / obj.frequency;
            timeHeader = strjoin(arrayfun(@(x) sprintf('%.4f', x), traceTime, 'UniformOutput', false), ', ');
            
            % Filter out out-of-range traces.
            nSamples = numel(obj.dff);
            k = triggerIds > halfWindow & triggerIds + halfWindow < nSamples;
            triggerLabels = triggerLabels(k);
            triggerIds = triggerIds(k);
            triggerData = obj.dff;
            triggerData = triggerData(window + triggerIds);
            triggerTime = obj.time(triggerIds);
            
            % All triggered windows with corresponding epoch label.
            % Order depends on epoch definitions. Overlapping is possible.
            % Rows represent a single trigger:
            % First column is the condition label of the trigger followed by the trace around each trigger, with each trigger at the center column (n / 2 + 1) labeled with "zero".
            fid = fopen(output1, 'w');
            fprintf(fid, ['# time, condition, ', timeHeader, '\n']);
            format = ['%.4f, %i', repmat(', %.4f', 1, triggeredWindow), '\n'];
            fprintf(fid, format, [triggerTime, triggerLabels, triggerData]');
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
            fprintf(fid, ['# condition, ', timeHeader, '\n']);
            format = ['%i', repmat(', %.4f', 1, triggeredWindow), '\n'];
            fprintf(fid, format, [uLabels, averages]');
            fclose(fid);
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

function plotTriggerAverage(data, ids, labels, window, frequency, colors, names)
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
        text(0.5, 0.5, 'No events', 'HorizontalAlignment', 'center');
    end
    legend('show');
    xlabel('Time (s)');
    axis('tight');
end