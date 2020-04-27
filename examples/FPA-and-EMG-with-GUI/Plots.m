% 2020-04-23. Leonardo Molina.
% 2020-04-23. Last modified.
classdef Plots < handle
    properties (Access = private)
        handles
        initialized
    end
    
    methods
        function obj = Plots()
            obj.initialized = false;
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
                    case 'threshold'
                        set(obj.handles.thresholdPlot, 'XData', time, 'YData', signal);
                    case 'xcorr'
                        set(obj.handles.xcorrPlot, 'XData', time, 'YData', signal);
                end
            end
        end
        
        function close(obj)
            delete(obj.handles.isolatedFigure);
            delete(obj.handles.mixedFigure);
        end
    end
    
    methods (Access = private)
        function setupFigures(obj)
            % Initialize and format plot handles.
            if ~obj.initialized || ~isobject(obj.handles.isolatedFigure)
                obj.handles.isolatedFigure = figure('name', 'FPA and EMG');
                obj.handles.fpaAxes = subplot(3, 1, 1);
                obj.handles.emgAxes = subplot(3, 1, 2);
                obj.handles.envelopeAxes = subplot(3, 1, 3);
                hold(obj.handles.emgAxes, 'all');
                
                obj.handles.fpaPlot = plot(obj.handles.fpaAxes, NaN(2, 1), NaN(2, 1), 'LineWidth', 1, 'DisplayName', 'Fiber-photometry');
                obj.handles.emgPlot = plot(obj.handles.emgAxes, NaN(2, 1), NaN(2, 1), 'DisplayName', 'Fiber-photometry');
                obj.handles.thresholdPlot = plot(obj.handles.emgAxes, NaN(2, 1), NaN(2, 1), 'LineWidth', 1, 'DisplayName', 'EMG threshold');
                obj.handles.envelopePlot1 = plot(obj.handles.emgAxes, NaN(2, 1), NaN(2, 1), 'LineWidth', 1, 'DisplayName', 'EMG envelope');
                obj.handles.envelopePlot2 = plot(obj.handles.envelopeAxes, NaN(2, 1), NaN(2, 1), 'LineWidth', 1, 'DisplayName', 'EMG envelope', 'Color', obj.handles.envelopePlot1.Color);
                linkaxes([obj.handles.fpaAxes, obj.handles.emgAxes, obj.handles.envelopeAxes], 'x');
                
                ylabel(obj.handles.fpaAxes, 'z-score');
                xlabel(obj.handles.emgAxes, 'time (s)');
                legend(obj.handles.fpaAxes, 'show');
                legend(obj.handles.emgAxes, 'show');
                title(obj.handles.fpaAxes, 'FP and EMG');
            end
            if ~obj.initialized || ~isobject(obj.handles.mixedFigure)
                obj.handles.mixedFigure = figure('name', 'FPA vs EMG');
                obj.handles.xcorrAxes = axes();
                obj.handles.xcorrPlot = plot(obj.handles.xcorrAxes, NaN(2, 1), NaN(2, 1), 'DisplayName', 'Cross-correlation');
                title(obj.handles.xcorrAxes, 'Cross-correlation between FP and EMG');
                xlabel(obj.handles.xcorrAxes, 'lag (s)');
                ylabel(obj.handles.xcorrAxes, 'xcorr');
            end
            
            obj.initialized = true;
        end
    end
end