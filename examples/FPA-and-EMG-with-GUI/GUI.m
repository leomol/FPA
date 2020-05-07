% 2020-04-23. Leonardo Molina.
% 2020-05-07. Last modified.
classdef GUI < handle
    properties (Dependent = true)
        envelopeSize
        envelopeLowpassFrequency
        envelopeThreshold
    end
    
    properties (Access = private)
        % Analysis settings.
        configuration
        
        % Results.
        emg
        fpa
        shared
        
        % Graphic handles.
        handles
        
        % Have figures ever been created.
        initialized
    end
    
    methods
        function obj = GUI(configuration)
            % Initialize.
            obj.configuration = configuration;
            obj.initialized = false;
            obj.emg = struct();
            obj.fpa = struct();
            obj.shared = struct();
            
            % Load data.
            obj.fpa.data = loadData(obj.configuration.fpa.file);
            obj.emg.data = loadABF(obj.configuration.emg.file);
            
            % Complete configuration.
            obj.configuration.fpa.bleachingEpochs = [configuration.conditionEpochs{2:2:end}];
            obj.configuration.fpa.resamplingFrequency = obj.configuration.resamplingFrequency;
            obj.configuration.fpa.conditionEpochs = obj.configuration.conditionEpochs;
            obj.configuration.fpa.plot = false;
            obj.configuration.emg.resamplingFrequency = obj.configuration.resamplingFrequency;
            obj.configuration.emg.conditionEpochs = obj.configuration.conditionEpochs;
            
            % Static analysis.
            obj.setupAnalysis();
            
            % Analysis.
            obj.recalculate();
            
            % Calculate valid ranges.
            obj.configuration.envelopeSizeRange = [0, 10];
            obj.configuration.envelopeLowpassFrequencyRange = [1e-6, obj.configuration.resamplingFrequency / 20];
            obj.configuration.envelopeThresholdRange = [0, prctile(obj.emg.signal, 99)];
            
            % Show figures and controls.
            obj.setupFigures();
            
            % Trigger updates.
            obj.envelopeSize = obj.envelopeSize;
            obj.envelopeLowpassFrequency = obj.envelopeLowpassFrequency;
            obj.envelopeThreshold = obj.envelopeThreshold;
            
            % Adjust axis limits once.
            axis([obj.handles.fpaAxes, obj.handles.emgAxes, obj.handles.envelopeAxes, obj.handles.xcorrAxes], 'tight');
        end
        
        function delete(obj)
            obj.close();
            delete(obj.handles.controlFigure);
        end
        
        function value = get.envelopeSize(obj)
            slider = obj.handles.envelopeSizeSlider;
            value = slider.Value;
        end
        
        function set.envelopeSize(obj, value)
            slider = obj.handles.envelopeSizeSlider;
            if value >= slider.Limits(1) && value <= slider.Limits(2)
                obj.configuration.emg.envelopeSize = value;
                slider.Value = value;
                obj.updateSlider(slider);
            else
                warn('.2f is outside the expected range (%.2f %.2f) for envelopeSize.', value, slider.Limits(1), slider.Limits(2));
            end
        end
        
        function value = get.envelopeLowpassFrequency(obj)
            slider = obj.handles.envelopeLowpassFrequencySlider;
            value = slider.Value;
        end
        
        function set.envelopeLowpassFrequency(obj, value)
            slider = obj.handles.envelopeLowpassFrequencySlider;
            if value >= slider.Limits(1) && value <= slider.Limits(2)
                obj.configuration.emg.envelopeLowpassFrequency = value;
                slider.Value = value;
                obj.updateSlider(slider);
            else
                warn('%.2f is outside the expected range (%.2f %.2f) for envelopeLowpassFrequency.', value, slider.Limits(1), slider.Limits(2));
            end
        end
        
        function value = get.envelopeThreshold(obj)
            slider = obj.handles.envelopeThresholdSlider;
            value = slider.Value;
        end
        
        function set.envelopeThreshold(obj, value)
            slider = obj.handles.envelopeThresholdSlider;
            if value >= slider.Limits(1) && value <= slider.Limits(2)
                obj.configuration.emg.envelopeThreshold = value;
                slider.Value = value;
                obj.updateSlider(slider);
            else
                warn('%.2f is outside the expected range (%.2f %.2f) for envelopeThreshold.', value, slider.Limits(1), slider.Limits(2));
            end
        end
    end
    
    methods (Access = private)
        function setupFigures(obj)
            % Initialize and format plot handles.
            if ~obj.initialized || ~isvalid(obj.handles.controlFigure)
                obj.handles.controlFigure = uifigure('name', sprintf('%s - control', mfilename('Class')), 'MenuBar', 'none', 'NumberTitle', 'off', 'ToolBar', 'none');
                obj.handles.controlFigure.CloseRequestFcn = @(handle, event)obj.close();

                dy = 50;
                nControls = 3;
                layout = uigridlayout(obj.handles.controlFigure);
                layout.RowHeight = repmat({dy}, [1, nControls]);
                layout.ColumnWidth = {'2x', '3x'};
                
                value = min(max(obj.configuration.emg.envelopeSize, obj.configuration.envelopeSizeRange(1)), obj.configuration.envelopeSizeRange(2));
                label = uilabel(layout);
                label.Layout.Column = 1;
                label.Layout.Row = 1;
                slider = uislider(layout, 'Limits', obj.configuration.envelopeSizeRange, 'Value', value, 'ValueChangedFcn', @(slider, event) obj.updateSlider(slider));
                slider.Layout.Column = 2;
                slider.Layout.Row = 1;
                obj.handles.envelopeSizeSlider = slider;
                obj.handles.envelopeSizeLabel = label;
                
                value = min(max(obj.configuration.emg.envelopeLowpassFrequency, obj.configuration.envelopeLowpassFrequencyRange(1)), obj.configuration.envelopeLowpassFrequencyRange(2));
                label = uilabel(layout);
                label.Layout.Column = 1;
                label.Layout.Row = 2;
                slider = uislider(layout, 'Limits', obj.configuration.envelopeLowpassFrequencyRange, 'Value', value, 'ValueChangedFcn', @(slider, event) obj.updateSlider(slider));
                slider.Layout.Column = 2;
                slider.Layout.Row = 2;
                obj.handles.envelopeLowpassFrequencySlider = slider;
                obj.handles.envelopeLowpassFrequencyLabel = label;
                
                value = min(max(obj.configuration.emg.envelopeThreshold, obj.configuration.envelopeThresholdRange(1)), obj.configuration.envelopeThresholdRange(2));
                label = uilabel(layout);
                label.Layout.Column = 1;
                label.Layout.Row = 3;
                slider = uislider(layout, 'Limits', obj.configuration.envelopeThresholdRange, 'Value', value, 'ValueChangedFcn', @(slider, event) obj.updateSlider(slider));
                slider.Layout.Column = 2;
                slider.Layout.Row = 3;
                obj.handles.envelopeThresholdSlider = slider;
                obj.handles.envelopeThresholdLabel = label;

                position = obj.handles.controlFigure.Position;
                obj.handles.controlFigure.Position = [position(1:3), dy * (nControls + 1)];
            end
            
            if ~obj.initialized || ~isvalid(obj.handles.rasterFigure)
                obj.handles.rasterFigure = figure('name', 'FPA and EMG');
                obj.handles.fpaAxes = subplot(3, 1, 1);
                obj.handles.emgAxes = subplot(3, 1, 2);
                obj.handles.envelopeAxes = subplot(3, 1, 3);
                hold(obj.handles.emgAxes, 'all');
                
                obj.handles.fpaPlot = plot(obj.handles.fpaAxes, NaN(2, 1), NaN(2, 1), 'DisplayName', 'Fiber-photometry');
                obj.handles.emgPlot = plot(obj.handles.emgAxes, NaN(2, 1), NaN(2, 1), 'DisplayName', 'EMG raw');
                obj.handles.envelopePlot1 = plot(obj.handles.emgAxes, NaN(2, 1), NaN(2, 1), 'LineWidth', 1.5, 'DisplayName', 'EMG envelope');
                obj.handles.envelopePlot2 = plot(obj.handles.envelopeAxes, NaN(2, 1), NaN(2, 1), 'DisplayName', 'EMG envelope', 'Color', obj.handles.envelopePlot1.Color);
                obj.handles.thresholdPlot = plot(obj.handles.emgAxes, NaN(2, 1), NaN(2, 1), 'LineStyle', '--', 'Color', [0, 0, 0], 'DisplayName', 'Envelope threshold');
                obj.handles.thresholdedPlot = plot(obj.handles.emgAxes, NaN(2, 1), NaN(2, 1), 'LineWidth', 1.5, 'DisplayName', 'Envelope thresholded');
                linkaxes([obj.handles.fpaAxes, obj.handles.emgAxes, obj.handles.envelopeAxes], 'x');
                
                ylabel(obj.handles.fpaAxes, 'z-score');
                xlabel(obj.handles.emgAxes, 'time (s)');
                legend(obj.handles.fpaAxes, 'show');
                legend(obj.handles.emgAxes, 'show');
                legend(obj.handles.envelopeAxes, 'show');
                title(obj.handles.fpaAxes, 'FP and EMG');
            end
            
            if ~obj.initialized || ~isvalid(obj.handles.xcorrFigure)
                obj.handles.xcorrFigure = figure('name', 'FPA vs EMG');
                obj.handles.xcorrAxes = axes();
                obj.handles.xcorrPlot = plot(obj.handles.xcorrAxes, NaN(2, 1), NaN(2, 1), 'DisplayName', 'Cross-correlation');
                title(obj.handles.xcorrAxes, 'Cross-correlation between FP and EMG');
                xlabel(obj.handles.xcorrAxes, 'lag (s)');
                ylabel(obj.handles.xcorrAxes, 'xcorr');
            end
            obj.initialized = true;
        end

        function updateSlider(obj, slider)
            switch slider
                case obj.handles.envelopeSizeSlider
                    baseText = 'Envelope size';
                    obj.configuration.emg.envelopeSize = slider.Value;
                    label = obj.handles.envelopeSizeLabel;
                case obj.handles.envelopeLowpassFrequencySlider
                    baseText = 'Envelope low-pass frequency';
                    obj.configuration.emg.envelopeLowpassFrequency = slider.Value;
                    label = obj.handles.envelopeLowpassFrequencyLabel;
                case obj.handles.envelopeThresholdSlider
                    baseText = 'Envelope threshold';
                    obj.configuration.emg.envelopeThreshold = slider.Value;
                    label = obj.handles.envelopeThresholdLabel;
            end
            label.Text = sprintf('%s (%.2f):', baseText, slider.Value);
            obj.recalculate();
            obj.refresh();
        end
        
        function recalculate(obj)
            % Envelope.
            envelopeSamples = max(2, ceil(obj.configuration.emg.envelopeSize * obj.configuration.resamplingFrequency));
            [high, low] = envelope(obj.emg.signal, envelopeSamples, 'rms');
            % Lowpass filter.
            lowpassFilter = designfilt('lowpassiir', 'HalfPowerFrequency', obj.configuration.emg.envelopeLowpassFrequency, 'SampleRate', obj.configuration.resamplingFrequency, 'DesignMethod', 'butter', 'FilterOrder', 12);
            high = filtfilt(lowpassFilter, high);
            low = filtfilt(lowpassFilter, low);
            obj.emg.envelope = [high, low];
            % Threshold envelope.
            obj.emg.thresholded = obj.configuration.emg.envelopeThreshold * (high >= obj.configuration.emg.envelopeThreshold);
            % Cross-correlation.
            xcLags = round(obj.configuration.xcorrSeconds * obj.configuration.resamplingFrequency);
            obj.shared.xcorrTics = (-xcLags:xcLags) / obj.configuration.resamplingFrequency;
            obj.shared.xcorr = xcorr(obj.fpa.signal, obj.emg.envelope(:, 1), xcLags, 'normalized');
        end
        
        function refresh(obj)
            obj.setData('fpa', obj.fpa.time, obj.fpa.signal, 'emg', obj.emg.time, obj.emg.signal, 'envelope', obj.emg.time, obj.emg.envelope(:, 1), 'thresholded', obj.emg.time, obj.emg.thresholded, 'threshold', obj.emg.time([1, end]), obj.configuration.emg.envelopeThreshold([1, 1]), 'xcorr', obj.shared.xcorrTics, obj.shared.xcorr);
        end
        
        function setData(obj, varargin)
            obj.setupFigures();
            targetNames = varargin(1:3:end);
            nTargets = numel(targetNames);
            for t = 1:nTargets
                targetName = targetNames{t};
                time = varargin{3 * t - 1};
                signal = varargin{3 * t - 0};
                switch targetName
                    case 'fpa'
                        set(obj.handles.fpaPlot, 'XData', time, 'YData', signal);
                    case 'emg'
                        set(obj.handles.emgPlot, 'XData', time, 'YData', signal);
                    case 'envelope'
                        set(obj.handles.envelopePlot1, 'XData', time, 'YData', signal);
                        set(obj.handles.envelopePlot2, 'XData', time, 'YData', signal);
                    case 'thresholded'
                        set(obj.handles.thresholdedPlot, 'XData', time, 'YData', signal);
                    case 'threshold'
                        set(obj.handles.thresholdPlot, 'XData', time, 'YData', signal);
                    case 'xcorr'
                        set(obj.handles.xcorrPlot, 'XData', time, 'YData', signal);
                end
            end
        end
        
        function setupAnalysis(obj)
            % Process FP.
            result = FPA(obj.fpa.data(:, 1), obj.fpa.data(:, obj.configuration.fpa.fp465Column), obj.fpa.data(:, obj.configuration.fpa.fp405Column), obj.configuration.fpa);
            fnames = fieldnames(result);
            nNames = numel(fnames);
            for i = 1:nNames
                fname = fnames{i};
                obj.fpa.(fname) = result.(fname);
            end
            
            obj.fpa.time = obj.fpa.time(obj.fpa.epochIds);
            obj.fpa.signal = obj.fpa.signal(obj.fpa.epochIds);
            obj.fpa.reference = obj.fpa.reference(obj.fpa.epochIds);
            obj.fpa.dff = obj.fpa.dff(obj.fpa.epochIds);
            % Print processing warnings.
            cellfun(@warning, obj.fpa.warnings);

            % Process EMG.
            obj.emg.warnings = {};
            % Band-pass (100 - 500 Hz), (rectify, resample, smooth ==> implemented with envelope function).
            obj.emg.time = obj.emg.data(:, 1);
            obj.emg.signal = obj.emg.data(:, obj.configuration.emg.emgColumn);
            obj.emg.sourceFrequency = 1 / mean(diff(obj.emg.time));
            % Highpass filter.
            bandpassFilter = designfilt('bandpassfir', 'CutoffFrequency1', obj.configuration.emg.bandpassFrequency(1), 'CutoffFrequency2', obj.configuration.emg.bandpassFrequency(2), 'SampleRate', obj.emg.sourceFrequency, 'DesignMethod', 'window', 'FilterOrder', 12);
            obj.emg.signal = filtfilt(bandpassFilter, obj.emg.signal);
            % Notch filter for 60 Hz.
            notchFilter = designfilt('bandstopiir', 'FilterOrder', 2, 'HalfPowerFrequency1', 59, 'HalfPowerFrequency2', 61, 'DesignMethod', 'butter', 'SampleRate', obj.emg.sourceFrequency);
            obj.emg.signal = filtfilt(notchFilter, obj.emg.signal);
            % Keep given time epochs.
            ids = time2id(obj.emg.time, [obj.configuration.conditionEpochs{2:2:end}]);
            obj.emg.time = obj.emg.time(ids);
            obj.emg.signal = obj.emg.signal(ids);
            % Detrend.
            obj.emg.signal = detrend(obj.emg.signal);
            % Resample to target frequency.
            if obj.configuration.resamplingFrequency < obj.emg.sourceFrequency
                % Express frequency as a ratio p/q.
                [p, q] = rat(obj.configuration.resamplingFrequency / obj.emg.sourceFrequency);
                % Resample: interpolate every p/q/f, upsample by p, filter, downsample by q.
                [obj.emg.signal, obj.emg.time] = resample(obj.emg.signal, obj.emg.time, obj.configuration.resamplingFrequency, p, q);
            else
                obj.emg.warnings{end + 1} = sprintf('[emg] Cannot resample to frequencies higher than the source frequency (%.2f Hz).', obj.emg.sourceFrequency);
            end

            % Print processing warnings.
            cellfun(@warning, obj.emg.warnings);
        end
        
        function close(obj)
            if obj.initialized && isvalid(obj.handles.rasterFigure)
                delete(obj.handles.rasterFigure);
            end
            if obj.initialized && isvalid(obj.handles.xcorrFigure)
                delete(obj.handles.xcorrFigure);
            end
            delete(obj.handles.controlFigure);
        end
    end
end

function warn(format, varargin)
    message = sprintf(format, varargin{:});
    warndlg(message, 'Warning', 'replace');
end