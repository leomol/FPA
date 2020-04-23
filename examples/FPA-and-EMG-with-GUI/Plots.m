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
                        target = obj.handles.fpaPlot;
                    case 'emg'
                        target = obj.handles.emgPlot;
                    case 'envelope'
                        target = obj.handles.envelopePlot;
                    case 'threshold'
                        target = obj.handles.thresholdPlot;
                    case 'xcorr'
                        target = obj.handles.xcorrPlot;
                end
                set(target, 'XData', time, 'YData', signal);
            end
            
        end
    end
    
    methods (Access = private)
        function setupFigures(obj)
            if ~obj.initialized || ~isobject(obj.handles.isolatedFigure)
                % Initialize and format plot handles.
                obj.handles.isolatedFigure = figure('name', 'FPA and EMG');
                obj.handles.fpaAxes = subplot(2, 1, 1);
                obj.handles.emgAxes = subplot(2, 1, 2);
                hold(obj.handles.emgAxes, 'all');

                obj.handles.fpaPlot = plot(obj.handles.fpaAxes, NaN(2, 1), NaN(2, 1), 'LineWidth', 1, 'DisplayName', 'Fiber-photometry');
                obj.handles.emgPlot = plot(obj.handles.emgAxes, NaN(2, 1), NaN(2, 1), 'DisplayName', 'EMG');
                obj.handles.envelopePlot = plot(obj.handles.emgAxes, NaN(2, 1), NaN(2, 1), 'LineWidth', 1, 'DisplayName', 'EMG envelope');
                obj.handles.thresholdPlot = plot(obj.handles.emgAxes, NaN(2, 1), NaN(2, 1), 'LineWidth', 1, 'DisplayName', 'EMG threshold');
                linkaxes([obj.handles.fpaAxes, obj.handles.emgAxes], 'x');

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