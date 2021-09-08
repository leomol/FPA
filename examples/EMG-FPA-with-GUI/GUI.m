% GUI for processing EMG and FP data.
% see FPA.

% 2020-11-01. Leonardo Molina.
% 2021-09-08. Last modified.
classdef GUI < handle
    properties
        fpa
    end
    
    properties % (Access = private)
        doricChannelsEntries
        behaviorEpochsEntries
        labChartChannelsEntries
        labChartBlocksEntries
        settingsFilename
        settings
        cache
        h
    end
    
    methods
        function obj = GUI(varargin)
            addpath('../../');
            addpath('../../common');
            
            obj.doricChannelsEntries = {};
            obj.behaviorEpochsEntries = {};
            obj.labChartChannelsEntries = {};
            obj.labChartBlocksEntries = {};
            obj.cache.doricFilename = '<empty>';
            obj.cache.labChartFilename = '<empty>';
            obj.cache.behaviorFilename = '<empty>';
            
            % Calculate GUI dimensions.
            screenSize = get(0, 'ScreenSize');
            screenSize = screenSize(3:4);
            targetSize = screenSize([1, 2]) .* [0.75, 0.75];
            uiHeight = 27;
            nRows = floor(targetSize(2) / uiHeight) - 2;
            
            % Create main figure.
            obj.h.control = uifigure('Name', 'EMG * FP analysis', 'MenuBar', 'none', 'NumberTitle', 'off', 'ToolBar', 'none', 'WindowState', 'maximized');
            obj.h.control.Position = [(screenSize(1) - targetSize(1)) / 2, (screenSize(2) - targetSize(2)) / 2, targetSize(1), targetSize(2)];
            
            % 3x3 layout.
            rowCount = floor(nRows / 3);
            mainLayout = uigridlayout(obj.h.control);
            mainLayout.RowHeight = {'1x', '1x', '0.5x'};
            mainLayout.ColumnWidth = repmat({'1x'}, 1, 3);
            
            % General panel.
            panel = uipanel(mainLayout, 'Title', 'General');
            panel = uigridlayout(panel, [1, 1]);
            panel.RowHeight = repmat({'1x'}, 1, rowCount);
            panel.ColumnWidth = {'3x', '2x'};
            
            % Resampling, behavior, xcorr.
            label = uilabel(panel, 'Text', 'Epochs type', 'HorizontalAlignment', 'right');
            label.Layout.Row = 1;
            label.Layout.Column = 1;
            label = uilabel(panel, 'Text', 'Resampling frequency', 'HorizontalAlignment', 'right');
            label.Layout.Row = 2;
            label.Layout.Column = 1;
            label = uilabel(panel, 'Text', 'Cross-correlation lag', 'HorizontalAlignment', 'right');
            label.Layout.Row = 3;
            label.Layout.Column = 1;
            
            % Epochs type.
            h = uidropdown(panel, 'Items', GUI.enumerate('epochsType'), 'ValueChangedFcn', @(h, ~)obj.onEpochsTypeDrop(h.Value));
            h.Layout.Row = 1;
            h.Layout.Column = 2;
            obj.h.epochsTypeDrop = h;
            h = uieditfield(panel, 'numeric', 'Limits', [10, Inf], 'ValueChangedFcn', @(h, ~)obj.saveSettings('resamplingFrequency', h.Value));
            h.Layout.Row = 2;
            h.Layout.Column = 2;
            obj.h.resamplingFrequencyEdit = h;
            h = uieditfield(panel, 'numeric', 'Limits', [1, Inf], 'ValueChangedFcn', @(h, ~)obj.saveSettings('xcorrLag', h.Value));
            h.Layout.Row = 3;
            h.Layout.Column = 2;
            obj.h.xcorrLagEdit = h;
            
            % Process button.
            h = uibutton(panel, 'Text', 'Process', 'ButtonPushedFcn', @(~, ~)obj.onProcess());
            h.Layout.Row = 5;
            h.Layout.Column = 2;
            obj.h.processButton = h;
            
            % EMG panel.
            panel = uipanel(mainLayout, 'Title', 'EMG');
            panel = uigridlayout(panel, [1, 1]);
            panel.RowHeight = repmat({'1x'}, 1, rowCount);
            panel.ColumnWidth = {'4x', '1x', '1x'};
            
            % EMG settings.
            label = uilabel(panel, 'Text', 'Bandpass type', 'HorizontalAlignment', 'right');
            label.Layout.Row = 1;
            label.Layout.Column = 1;
            label = uilabel(panel, 'Text', 'Envelope type', 'HorizontalAlignment', 'right');
            label.Layout.Row = 2;
            label.Layout.Column = 1;
            label = uilabel(panel, 'Text', 'Envelope lowpass frequency (Hz)', 'HorizontalAlignment', 'right');
            label.Layout.Row = 3;
            label.Layout.Column = 1;
            label = uilabel(panel, 'Text', 'Envelope size', 'HorizontalAlignment', 'right');
            label.Layout.Row = 4;
            label.Layout.Column = 1;
            label = uilabel(panel, 'Text', 'Bandpass frequency', 'HorizontalAlignment', 'right');
            label.Layout.Row = 5;
            label.Layout.Column = 1;
            
            h = uidropdown(panel, 'Items', GUI.enumerate('emgBandpassType'), 'ValueChangedFcn', @(h, ~)obj.saveSettings('emgBandpassType', h.Value));
            h.Layout.Row = 1;
            h.Layout.Column = [2, 3];
            obj.h.emgBandpassTypeDrop = h;
            h = uidropdown(panel, 'Items', GUI.enumerate('emgEnvelopeType'), 'ValueChangedFcn', @(h, ~)obj.onEmgEnvelopeTypeDrop(h.Value));
            h.Layout.Row = 2;
            h.Layout.Column = [2, 3];
            obj.h.emgEnvelopeTypeDrop = h;
            h = uieditfield(panel, 'numeric', 'Limits', [1, Inf], 'ValueChangedFcn', @(h, ~)obj.saveSettings('emgEnvelopeLowpassFrequency', h.Value));
            h.Layout.Row = 3;
            h.Layout.Column = [2, 3];
            obj.h.emgEnvelopeLowpassFrequencyEdit = h;
            h = uieditfield(panel, 'numeric', 'Limits', [0.1, Inf], 'ValueChangedFcn', @(h, ~)obj.saveSettings('emgEnvelopeSize', h.Value));
            h.Layout.Row = 4;
            h.Layout.Column = [2, 3];
            obj.h.emgEnvelopeSizeEdit = h;
            h = uieditfield(panel, 'numeric', 'Limits', [0.1, Inf], 'ValueChangedFcn', @(h, ~)obj.saveSettings('emgBandpassFrequencyLow', h.Value));
            h.Layout.Row = 5;
            h.Layout.Column = 2;
            obj.h.emgBandpassFrequencyLowEdit = h;
            h = uieditfield(panel, 'numeric', 'Limits', [0.1, Inf], 'ValueChangedFcn', @(h, ~)obj.saveSettings('emgBandpassFrequencyHigh', h.Value));
            h.Layout.Row = 5;
            h.Layout.Column = 3;
            obj.h.emgBandpassFrequencyHighEdit = h;
            
            % FP panel.
            panel = uipanel(mainLayout, 'Title', 'Fiber-photometry');
            panel = uigridlayout(panel, [1, 1]);
            panel.RowHeight = repmat({'1x'}, 1, rowCount);
            panel.ColumnWidth = {'4x', '2x'};
            
            % FP settings.
            label = uilabel(panel, 'Text', 'Baseline correction type', 'HorizontalAlignment', 'right');
            label.Layout.Row = 1;
            label.Layout.Column = 1;
            label = uilabel(panel, 'Text', 'Baseline lowpass frequency', 'HorizontalAlignment', 'right');
            label.Layout.Row = 2;
            label.Layout.Column = 1;
            label = uilabel(panel, 'Text', 'Baseline epochs', 'HorizontalAlignment', 'right');
            label.Layout.Row = 3;
            label.Layout.Column = 1;
            label = uilabel(panel, 'Text', 'Fit reference', 'HorizontalAlignment', 'right');
            label.Layout.Row = 4;
            label.Layout.Column = 1;
            
            label = uilabel(panel, 'Text', 'Artifact epochs', 'HorizontalAlignment', 'right');
            label.Layout.Row = 5;
            label.Layout.Column = 1;
            label = uilabel(panel, 'Text', 'Peaks lowpass frequency', 'HorizontalAlignment', 'right');
            label.Layout.Row = 6;
            label.Layout.Column = 1;
            label = uilabel(panel, 'Text', 'df/f lowpass frequency', 'HorizontalAlignment', 'right');
            label.Layout.Row = 7;
            label.Layout.Column = 1;
            label = uilabel(panel, 'Text', 'Thresholding function', 'HorizontalAlignment', 'right');
            label.Layout.Row = 8;
            label.Layout.Column = 1;
            label = uilabel(panel, 'Text', 'Thresholding factor', 'HorizontalAlignment', 'right');
            label.Layout.Row = 9;
            label.Layout.Column = 1;
            label = uilabel(panel, 'Text', 'Normalization type', 'HorizontalAlignment', 'right');
            label.Layout.Row = 10;
            label.Layout.Column = 1;
            
            h = uidropdown(panel, 'Items', GUI.enumerate('fpBaselineType'), 'ValueChangedFcn', @(h, ~)obj.onFPBaselineTypeDrop(h.Value));
            h.Layout.Row = 1;
            h.Layout.Column = 2;
            obj.h.fpBaselineTypeDrop = h;
            h = uieditfield(panel, 'numeric', 'Limits', [1e-3, Inf], 'ValueChangedFcn', @(h, ~)obj.saveSettings('fpBaselineLowpassFrequency', h.Value));
            h.Layout.Row = 2;
            h.Layout.Column = 2;
            obj.h.fpBaselineLowpassFrequencyEdit = h;
            h = uieditfield(panel, 'ValueChangedFcn', @(~, ~)obj.onBaselineEpochsEdit);
            h.Layout.Row = 3;
            h.Layout.Column = 2;
            obj.h.fpBaselineEpochsEdit = h;
            h = uidropdown(panel, 'Items', GUI.enumerate('fpFitReference'), 'ValueChangedFcn', @(h, ~)obj.saveSettings('fpFitReference', h.Value));
            h.Layout.Row = 4;
            h.Layout.Column = 2;
            obj.h.fpFitReferenceDrop = h;
            
            h = uieditfield(panel, 'ValueChangedFcn', @(~, ~)obj.onArtifactEpochsEdit);
            h.Layout.Row = 5;
            h.Layout.Column = 2;
            obj.h.fpArtifactEpochsEdit = h;
            h = uieditfield(panel, 'numeric', 'Limits', [1e-3, Inf], 'ValueChangedFcn', @(h, ~)obj.saveSettings('fpPeaksLowpassFrequency', h.Value));
            h.Layout.Row = 6;
            h.Layout.Column = 2;
            obj.h.fpPeaksLowpassFrequencyEdit = h;
            h = uieditfield(panel, 'numeric', 'Limits', [1e-3, Inf], 'ValueChangedFcn', @(h, ~)obj.saveSettings('fpLowpassFrequency', h.Value));
            h.Layout.Row = 7;
            h.Layout.Column = 2;
            obj.h.fpLowpassFrequencyEdit = h;
            h = uidropdown(panel, 'Items', GUI.enumerate('fpThresholdingFunctions'), 'ValueChangedFcn', @(h, ~)obj.saveSettings('fpThresholdingFunctions', h.Value));
            h.Layout.Row = 8;
            h.Layout.Column = 2;
            obj.h.fpThresholdingFunctionsDrop = h;
            h = uieditfield(panel, 'numeric', 'Limits', [1e-3, Inf], 'ValueChangedFcn', @(h, ~)obj.saveSettings('fpThresholdingFactor', h.Value));
            h.Layout.Row = 9;
            h.Layout.Column = 2;
            obj.h.fpThresholdingFactorEdit = h;
            h = uidropdown(panel, 'Items', GUI.enumerate('fpNormalizationType'), 'ValueChangedFcn', @(h, ~)obj.saveSettings('fpNormalizationType', h.Value));
            h.Layout.Row = 10;
            h.Layout.Column = 2;
            obj.h.fpNormalizationTypeDrop = h;
            
            % Epochs panel.
            panel = uipanel(mainLayout, 'Title', 'Manual Epoch definitions');
            obj.h.conditionEpochsPanel = panel;
            panel = uigridlayout(panel, [1, 1]);
            panel.RowHeight = {'1x'};
            panel.ColumnWidth = {'1x'};
            h = uitextarea(panel, 'ValueChangedFcn', @(~, ~)obj.onConditionEpochsEdit);
            obj.h.conditionEpochsEdit = h;
            
            % CleverSys epochs panel.
            panel = uipanel(mainLayout, 'Title', 'Epoch definitions');
            panel.Layout = obj.h.conditionEpochsPanel.Layout;
            obj.h.behaviorPanel = panel;
            panel = uigridlayout(panel, [1, 1]);
            panel.RowHeight = repmat({'1x'}, 1, rowCount);
            panel.ColumnWidth = {'5x', '1x'};
            
            % CleverSys filename.
            h = uieditfield(panel, 'ValueChangedFcn', @(~, ~)obj.onCleverSysFilenameEdit);
            h.Layout.Row = 1;
            h.Layout.Column = 1;
            obj.h.behaviorFilenameEdit = h;
            h = uibutton(panel, 'Text', 'Select', 'ButtonPushedFcn', @(~, ~)obj.onGetFile(obj.h.behaviorFilenameEdit));
            h.Layout.Row = 1;
            h.Layout.Column = 2;
            obj.h.behaviorButton = h;
            
            % CleverSys sheet dropdown and list.
            h = uidropdown(panel, 'Items', {}, 'ValueChangedFcn', @(h, ~)obj.onCleverSysSheetDrop(h.Value));
            h.Layout.Row = 2;
            h.Layout.Column = [1, 2];
            obj.h.cleverSysSheetDrop = h;
            h = uilistbox(panel, 'Items', {}, 'Multiselect', true, 'FontName', 'Monospaced');
            h.Layout.Row = [3, rowCount];
            h.Layout.Column = [1, 2];
            obj.addMenu(h, {'Use', 'Use', 'Unset', 'Unset'}, @obj.setChoice);
            obj.h.behaviorEpochsList = h;
            
            % Doric panel.
            panel = uipanel(mainLayout, 'Title', 'Doric Channels');
            panel = uigridlayout(panel, [1, 1]);
            panel.RowHeight = repmat({'1x'}, 1, rowCount);
            panel.ColumnWidth = {'5x', '1x'};
            
            % Doric filename edit, button and channel list.
            h = uieditfield(panel, 'ValueChangedFcn', @(~, ~)obj.onDoricFilenameEdit);
            h.Layout.Row = 1;
            h.Layout.Column = 1;
            obj.h.doricFilenameEdit = h;
            h = uibutton(panel, 'Text', 'Select', 'ButtonPushedFcn', @(~, ~)obj.onGetFile(obj.h.doricFilenameEdit));
            h.Layout.Row = 1;
            h.Layout.Column = 2;
            obj.h.doricButton = h;
            
            h = uilistbox(panel, 'Items', {}, 'Multiselect', false, 'FontName', 'Monospaced');
            h.Layout.Row = [2, rowCount];
            h.Layout.Column = [1, 2];
            obj.addMenu(h, {'Set as time', 'Time', 'Set as signal', 'Signal', 'Set as reference', 'Reference', 'Set as camera', 'Camera', 'Unset', 'Unset'}, @obj.setChoice);
            obj.h.doricChannelsList = h;
            
            % LabChart filename edit, button and channel list.
            panel = uipanel(mainLayout, 'Title', 'LabChart channels');
            panel = uigridlayout(panel, [1, 1]);
            panel.RowHeight = repmat({'1x'}, 1, rowCount);
            panel.ColumnWidth = {'5x', '1x'};
            h = uieditfield(panel, 'ValueChangedFcn', @(~, ~)obj.onLabChartFilenameEdit);
            h.Layout.Row = 1;
            h.Layout.Column = 1;
            obj.h.labChartFilenameEdit = h;
            h = uibutton(panel, 'Text', 'Select', 'ButtonPushedFcn', @(~, ~)obj.onGetFile(obj.h.labChartFilenameEdit));
            h.Layout.Row = 1;
            h.Layout.Column = 2;
            obj.h.labChartButton = h;
            h = uilistbox(panel, 'Items', {}, 'Multiselect', false, 'FontName', 'Monospaced');
            h.Layout.Row = [2, rowCount];
            h.Layout.Column = [1, 2];
            obj.addMenu(h, {'Set as EMG', 'EMG', 'Unset', 'Unset'}, @obj.setChoice);
            obj.h.labChartChannelsList = h;
            
            % LabChart/Blocks panel.
            panel = uipanel(mainLayout, 'Title', 'LabChart blocks');
            panel = uigridlayout(panel, [1, 1]);
            panel.RowHeight = repmat({'1x'}, 1, rowCount);
            panel.ColumnWidth = {'1x'};
            h = uilistbox(panel, 'Items', {}, 'Multiselect', false, 'FontName', 'Monospaced');
            h.Layout.Row = [1, rowCount];
            h.Layout.Column = 1;
            obj.addMenu(h, {'Use this block', 'Use', 'Unset', 'Unset'}, @obj.setChoice);
            obj.h.labChartBlocksList = h;
            
            % LabChart/Comments panel.
            panel = uipanel(mainLayout, 'Title', 'LabChart comments');
            panel = uigridlayout(panel, [1, 1]);
            panel.RowHeight = repmat({'1x'}, 1, rowCount);
            panel.ColumnWidth = {'1x'};
            h = uihtml(panel);
            h.Layout.Row = [1, rowCount];
            h.Layout.Column = 1;
            h.HTMLSource = '<div></div>';
            obj.h.labChartCommentsHtml = h;
            
            % Plots panel.
            panel = uipanel(mainLayout, 'Title', 'Plots');
            panel = uigridlayout(panel, [1, 1]);
            panel.RowHeight = repmat({'1x'}, 1, ceil(rowCount / 2));
            panel.ColumnWidth = {'1x', '1x'};
            
            savePlotSelection = @()obj.saveSettings('plotFpTrace', obj.h.fpPlots(1).Value, 'plotFpPower', obj.h.fpPlots(2).Value, 'plotFpStats', obj.h.fpPlots(3).Value, 'plotFpTrigger', obj.h.fpPlots(4).Value, 'plotFpTriggerAverage', obj.h.fpPlots(5).Value, 'plotFpAUC', obj.h.fpPlots(6).Value, 'plotEmgTrace', obj.h.emgPlots(1).Value, 'plotXcorr', obj.h.generalPlots(1).Value);
            obj.h.fpPlots(1) = uicheckbox(panel, 'Text', 'df/f trace and peaks', 'ValueChangedFcn', @(~, ~)savePlotSelection());
            obj.h.fpPlots(2) = uicheckbox(panel, 'Text', 'df/f power spectrum', 'ValueChangedFcn', @(~, ~)savePlotSelection());
            obj.h.fpPlots(3) = uicheckbox(panel, 'Text', 'df/f stats boxplot', 'ValueChangedFcn', @(~, ~)savePlotSelection());
            obj.h.fpPlots(4) = uicheckbox(panel, 'Text', 'df/f triggered', 'ValueChangedFcn', @(~, ~)savePlotSelection());
            obj.h.fpPlots(5) = uicheckbox(panel, 'Text', 'df/f triggered average', 'ValueChangedFcn', @(~, ~)savePlotSelection());
            obj.h.fpPlots(6) = uicheckbox(panel, 'Text', 'df/f normalized AUC', 'ValueChangedFcn', @(~, ~)savePlotSelection());
            obj.h.emgPlots(1) = uicheckbox(panel, 'Text', 'EMG trace and envelope', 'ValueChangedFcn', @(~, ~)savePlotSelection());
            obj.h.generalPlots(1) = uicheckbox(panel, 'Text', 'Cross-correlation df/f to EMG', 'ValueChangedFcn', @(~, ~)savePlotSelection());
            
            if nargin == 1
                switch class(varargin{1})
                    case {'char', 'string'}
                        obj.settingsFilename = varargin{1};
                        settings = load(obj.settingsFilename);
                        settings = settings.settings;
                        fprintf('Loaded settings from "%s"\n', obj.settingsFilename);
                    case 'struct'
                        settings = varargin{1};
                end
            elseif nargin == 0
                obj.settingsFilename = fullfile(getHome(), 'Documents', 'MATLAB', 'csmopto-settings.mat');
                if exist(obj.settingsFilename, 'file') == 2
                    settings = load(obj.settingsFilename);
                    settings = settings.settings;
                    fprintf('Loaded settings from "%s"\n', obj.settingsFilename);
                else
                    settings = getDefaults();
                end
            end
            defaults = getDefaults();
            obj.applySettings(defaults, settings);
        end
    end
    
    methods % (Access = private)
        function saveSettings(obj, varargin)
            keys = varargin(1:2:end);
            values = varargin(2:2:end);
            nKeys = numel(keys);
            settings = obj.settings;
            for i = 1:nKeys
                key = keys{i};
                value = values{i};
                settings.(key) = value;
            end
            obj.settings = settings;
            save(obj.settingsFilename, 'settings');
        end
        
        function applySettings(obj, defaults, settings)
            expectedNames = fieldnames(defaults);
            for i = 1:numel(expectedNames)
                name = expectedNames{i};
                if isfield(settings, name)
                    options = GUI.enumerate(name);
                    missingValue = false;
                    if isempty(options) || ismember(settings.(name), options)
                        wrongValue = false;
                    else
                        wrongValue = true;
                    end
                else
                    missingValue = true;
                    wrongValue = false;
                end
                if missingValue
                    value = defaults.(name);
                    warning('[GUI] "%s" is missing, using default:\n%s\n\n', name, strtrim(evalc('disp(value)')));
                elseif wrongValue
                    value = defaults.(name);
                    warning('[GUI] Provided value for "%s" is incorrect, using default:\n%s\n\n', name, strtrim(evalc('disp(value)')));
                else
                    value = settings.(name);
                end
                obj.settings.(name) = value;
            end
            
            if exist(settings.lastFolder, 'dir') ~= 7
                warning('[GUI] lastFolder "%s" does not exist, using default: "%s"', settings.lastFolder, defaults.lastFolder);
                obj.settings.lastFolder = defaults.lastFolder;
            end
            
            unexpectedNames = setdiff(fieldnames(settings), expectedNames);
            if numel(unexpectedNames) > 0
                warning('[GUI] Unexpected settings in "%s":%s', obj.settingsFilename, sprintf(' "%s"', unexpectedNames{:}));
            end
            
            obj.h.resamplingFrequencyEdit.Value = obj.settings.resamplingFrequency;
            obj.h.xcorrLagEdit.Value = obj.settings.xcorrLag;
            
            obj.h.emgBandpassTypeDrop.Value = obj.settings.emgBandpassType;
            obj.h.emgEnvelopeLowpassFrequencyEdit.Value = obj.settings.emgEnvelopeLowpassFrequency;
            obj.h.emgEnvelopeSizeEdit.Value = obj.settings.emgEnvelopeSize;
            obj.h.emgBandpassFrequencyLowEdit.Value = obj.settings.emgBandpassFrequencyLow;
            obj.h.emgBandpassFrequencyHighEdit.Value = obj.settings.emgBandpassFrequencyHigh;
            
            obj.h.fpBaselineTypeDrop.Value = obj.settings.fpBaselineType;
            obj.h.fpBaselineLowpassFrequencyEdit.Value = obj.settings.fpBaselineLowpassFrequency;
            obj.h.fpBaselineEpochsEdit.Value = obj.settings.fpBaselineEpochsText;
            obj.h.fpFitReferenceDrop.Value = obj.settings.fpFitReference;
            
            obj.h.fpPeaksLowpassFrequencyLowEdit.Value = obj.settings.fpPeaksLowpassFrequency;
            obj.h.fpLowpassFrequencyEdit.Value = obj.settings.fpLowpassFrequency;
            obj.h.fpThresholdingFunctionsDrop.Value = obj.settings.fpThresholdingFunctions;
            obj.h.fpThresholdingFactorEdit.Value = obj.settings.fpThresholdingFactor;
            obj.h.fpNormalizationTypeDrop.Value = obj.settings.fpNormalizationType;
            
            obj.h.conditionEpochsEdit.Value = obj.settings.conditionEpochsText;
            obj.h.fpArtifactEpochsEdit.Value = obj.settings.fpArtifactEpochsText;
            
            obj.h.behaviorFilenameEdit.Value = obj.settings.behaviorFilename;
            obj.h.doricFilenameEdit.Value = obj.settings.doricFilename;
            obj.h.labChartFilenameEdit.Value = obj.settings.labChartFilename;
            
            obj.h.fpPlots(1).Value = obj.settings.plotFpTrace;
            obj.h.fpPlots(2).Value = obj.settings.plotFpPower;
            obj.h.fpPlots(3).Value = obj.settings.plotFpStats;
            obj.h.fpPlots(4).Value = obj.settings.plotFpTrigger;
            obj.h.fpPlots(5).Value = obj.settings.plotFpTriggerAverage;
            obj.h.fpPlots(6).Value = obj.settings.plotFpAUC;
            obj.h.emgPlots(1).Value = obj.settings.plotEmgTrace;
            obj.h.generalPlots(1).Value = obj.settings.plotXcorr;
            
            % Cascade calls.
            obj.onDoricFilenameEdit();
            obj.onCleverSysFilenameEdit();
            obj.onLabChartFilenameEdit();
            obj.onFPBaselineTypeDrop(obj.settings.fpBaselineType);
            obj.onEpochsTypeDrop(obj.settings.epochsType);
            obj.onEmgEnvelopeTypeDrop(obj.settings.emgEnvelopeType);
        end
        
        function onProcess(obj)
            % Large data files exported from LabChart require an SDK for loading into MATLAB (it is not the typical memory limitation).
            % Using Doric's demodulated-FP
            % Ignoring Powerlab's raw-FP data.
            % Behavior is aligned with FP and EMG using the 20 TTL pulses occurring at 600s.
            % Change recording setup to save camera frame triggers.
            % 
            % FP1/FP2/EMG: continuosly from the start.
            % Camera trigger: Once at the start.
            % TTL sync input: 10min after start, Square wave (20 pulses, for 9.5s)
            
            set(obj.h.processButton, 'Text', 'Processing...');
            % !! set(obj.h.processButton, 'Enable', false);
            drawnow();
            messages = {};
            emgRequired = obj.settings.labChartFilename ~= "";
            
            % Fiber-photometry settings.
            fpTimeChannel = getIndex(obj.doricChannelsEntries, obj.settings.doricChannelsChoices, 'Time');
            fpSignalChannel = getIndex(obj.doricChannelsEntries, obj.settings.doricChannelsChoices, 'Signal');
            fpReferenceChannel = getIndex(obj.doricChannelsEntries, obj.settings.doricChannelsChoices, 'Reference');
            fpCameraChannel = getIndex(obj.doricChannelsEntries, obj.settings.doricChannelsChoices, 'Camera');
            if numel(fpTimeChannel) ~= 1
                messages = cat(2, messages, 'Select a single Time channel in Doric Channels.');
            end
            if numel(fpSignalChannel) ~= 1
                messages = cat(2, messages, 'Select a single Signal channel in Doric Channels.');
            end
            if numel(fpReferenceChannel) ~= 1
                messages = cat(2, messages, 'Select a single Reference channel in Doric Channels.');
            end
            if emgRequired && numel(fpCameraChannel) ~= 1
                messages = cat(2, messages, 'Select a single Camera channel in Doric Channels.');
            end
            if numel(messages) > 0
                errorDialog(messages);
                set(obj.h.processButton, 'Enable', true, 'Text', 'Process');
                return;
            end
            
            fp = struct();
            fp.baselineEpochs = obj.settings.fpBaselineEpochs;
            if isempty(fp.baselineEpochs)
                fp.baselineEpochs = [-Inf, Inf];
            end
            fp.artifactEpochs = obj.settings.fpArtifactEpochs;
            if isempty(fp.artifactEpochs)
                fp.artifactEpochs = [];
            end
            
            fp.airPLS = obj.settings.fpBaselineType == "airPLS";
            fp.fitReference = obj.settings.fpFitReference == "true";
            fp.baselineLowpassFrequency = obj.settings.fpBaselineLowpassFrequency;
            fp.peaksLowpassFrequency = obj.settings.fpPeaksLowpassFrequency;
            fp.lowpassFrequency = obj.settings.fpLowpassFrequency;
            fp.resamplingFrequency = obj.settings.resamplingFrequency;
            switch obj.settings.fpThresholdingFunctions
                case 'median:mad'
                    fp.threshold = {obj.settings.fpThresholdingFactor, @mad, @median};
                case 'mean:std'
                    fp.threshold = {obj.settings.fpThresholdingFactor, @std, @mean};
            end
            
            switch obj.settings.fpNormalizationType
                case 'z-score'
                    fp.f0 = @mean;
                    fp.f1 = @std;
                case 'df/f'
                    fp.f0 = @mean;
                    fp.f1 = @mean;
                case 'median:mad'
                    fp.f0 = @median;
                    fp.f1 = @mad;
                case 'median:std'
                    fp.f0 = @median;
                    fp.f1 = @std;
                case 'mean:mad'
                    fp.f0 = @mean;
                    fp.f1 = @mad;
            end
            
            % LabChart settings.
            if emgRequired
                labChart = struct();
                labChart.emgChannel = getIndex(obj.labChartChannelsEntries, obj.settings.labChartChannelsChoices, 'EMG');
                labChart.block = getIndex(obj.labChartBlocksEntries, obj.settings.labChartBlocksChoices, 'Use');
                if numel(labChart.emgChannel) ~= 1
                    messages = cat(2, messages, 'Select a single EMG channel in LabChart Channels.');
                end
                if numel(labChart.block) ~= 1
                    messages = cat(2, messages, 'Select a single Block in LabChart Blocks.');
                end
                if numel(messages) > 0
                    errorDialog(messages);
                    set(obj.h.processButton, 'Enable', true, 'Text', 'Process');
                    return;
                end
            end
            
            fprintf('Loading Doric data ... ');
            if isequal(obj.settings.doricFilename, obj.cache.doricFilename)
                fprintf('reused cache.\n');
                fpData = obj.cache.fpData;
            else
                fprintf('\n');
                fpData = loadDoric(obj.settings.doricFilename);
                obj.cache.fpData = fpData;
                obj.cache.doricFilename = obj.settings.doricFilename;
            end
            
            fpTime = fpData(:, fpTimeChannel);
            fpSignal = fpData(:, fpSignalChannel);
            fpReference = fpData(:, fpReferenceChannel);
            
            if emgRequired
                fprintf('Loading LabChart data ... ');
                if isequal(obj.settings.labChartFilename, obj.cache.labChartFilename)
                    fprintf('reused cache.\n');
                    labChartTime = obj.cache.labChartTime;
                    labChartData = obj.cache.labChartData;
                else
                    fprintf('\n');
                    [labChartTime, labChartData] = Aditch.getData(obj.settings.labChartFilename);
                    obj.cache.labChartTime = labChartTime;
                    obj.cache.labChartData = labChartData;
                    obj.cache.labChartFilename = obj.settings.labChartFilename;
                end
                emgTime = labChartTime{labChart.emgChannel, labChart.block};
                emgSignal = labChartData{labChart.emgChannel, labChart.block};
            end
            
            fprintf('Processing data ...\n');
            if emgRequired
                [mn1, mx1] = bounds(fpTime);
                [mn2, mx2] = bounds(emgTime);
                mn = max(mn1, mn2);
                mx = min(mx1, mx2);
                k1 = fpTime >= mn & fpTime <= mx;
                fpTime = fpTime(k1);
                fpSignal = fpSignal(k1);
                fpReference = fpReference(k1);
                k2 = emgTime >= mn & emgTime <= mx;
                emgTime = emgTime(k2);
                emgSignal = emgSignal(k2);
            
                [emgTime, emgSignal] = alignTimestamps(emgTime, 1 / obj.settings.resamplingFrequency, emgSignal);
                emgSourceFrequency = 1 / median(diff(emgTime));
                
                % Bandpass emg to remove artifacts.
                switch obj.settings.emgBandpassType
                    case 'zero-phase'
                        % High-order, butter filter, without phase shift.
                        bandpassFilter = designfilt('bandpassiir', 'HalfPowerFrequency1', obj.settings.emgBandpassFrequency(1), 'HalfPowerFrequency2', obj.settings.emgBandpassFrequency(2), 'SampleRate', emgSourceFrequency, 'DesignMethod', 'butter', 'FilterOrder', 12);
                        emgSignal = filtfilt(bandpassFilter, emgSignal);
                    case 'butter'
                        % Standard butter filter.
                        % http://www1.udel.edu/biology/rosewc/kaap686/notes/EMG%20analysis.pdf
                        fn = emgSourceFrequency / 2;
                        fl = obj.settings.emgBandpassFrequency(1) / fn;
                        fh = obj.settings.emgBandpassFrequency(2) / fn;
                        [b, a] = butter(4, [fl, fh], 'bandpass');
                        emgSignal = filter(b, a, emgSignal);
                    case 'window'
                        % MATLAB's default filter.
                        bandpassFilter = designfilt('bandpassfir', 'CutoffFrequency1', obj.settings.emgBandpassFrequency(1), 'CutoffFrequency2', obj.settings.emgBandpassFrequency(2), 'SampleRate', emgSourceFrequency, 'DesignMethod', 'window', 'FilterOrder', 12);
                        emgSignal = filtfilt(bandpassFilter, emgSignal);
                end

                % Notch filter to remove 60Hz noise.
                notchFilter = designfilt('bandstopiir', 'FilterOrder', 2, 'HalfPowerFrequency1', 59, 'HalfPowerFrequency2', 61, 'DesignMethod', 'butter', 'SampleRate', emgSourceFrequency);
                emgSignal = filtfilt(notchFilter, emgSignal);
                % Detrend.
                emgSignal = detrend(emgSignal);
                % Make data have the same sampling rate.
                if obj.settings.resamplingFrequency < emgSourceFrequency
                    [p, q] = rat(obj.settings.resamplingFrequency / emgSourceFrequency);
                    [emgSignal, emgTime] = resample(emgSignal, emgTime, obj.settings.resamplingFrequency, p, q);
                else
                    warning('[emg] Cannot resample to frequencies higher than the source frequency (%.2f Hz).', emgSourceFrequency);
                end
                switch obj.settings.emgEnvelopeType
                    case 'Low-pass'
                        % Rectify then 100Hz low pass filter.
                        lowpassFilter = designfilt('lowpassiir', 'HalfPowerFrequency', obj.settings.emgEnvelopeLowpassFrequency, 'SampleRate', obj.settings.resamplingFrequency, 'DesignMethod', 'butter', 'FilterOrder', 12);
                        emgHigh = filtfilt(lowpassFilter, abs(emgSignal));
                    case 'RMS'
                        % RMS.
                        envelopeSamples = max(2, ceil(obj.settings.emgEnvelopeSize * obj.settings.resamplingFrequency));
                        emgHigh = envelope(abs(emgSignal), envelopeSamples, 'rms');
                end
            end
            
            manualEpochs = obj.settings.conditionEpochs;
            if ~isempty(obj.settings.behaviorFilename) && exist(obj.settings.behaviorFilename, 'file') == 2
                fprintf('Loading epochs data ... ');
                if isequal(obj.settings.behaviorFilename, obj.cache.behaviorFilename) && isequal(obj.settings.cleverSysSheet, obj.cache.cleverSysSheet)
                    fprintf('reused cache.\n');
                    behaviorEpochs = obj.cache.behaviorEpochs;
                else
                    fprintf('\n');
                    behaviorEpochs = loadBehavior(obj.settings.behaviorFilename, obj.settings.cleverSysSheet);
                    obj.cache.behaviorEpochs = behaviorEpochs;
                    obj.cache.behaviorFilename = obj.settings.behaviorFilename;
                    obj.cache.cleverSysSheet = obj.settings.cleverSysSheet;
                end
                k = getIndex(obj.behaviorEpochsEntries, obj.settings.behaviorEpochsChoices, 'Use');
                k = sort([2 * k, 2 * k - 1]);
                behaviorEpochs = behaviorEpochs(k);
                epochs = [manualEpochs, behaviorEpochs];
            else
                epochs = manualEpochs;
            end
            if isempty(epochs)
                epochs = {'everything', [-Inf, Inf]};
            end
            
            % FP and EMG alignment.
            if emgRequired
                fpCameraData = fpData(:, fpCameraChannel);
                behaviorStart = fpTime(find(fpCameraData, 1));
                epochs(2:2:end) = cellfun(@(epoch) epoch + behaviorStart, epochs(2:2:end), 'UniformOutput', false);
            end
            
            fp.conditionEpochs = epochs;
            [fpTime, fpSignal, fpReference] = alignTimestamps(fpTime, 1 / obj.settings.resamplingFrequency, fpSignal, fpReference);
            obj.fpa = FPA(fpTime, fpSignal, fpReference, fp);
            cellfun(@warning, obj.fpa.warnings);
            
            if obj.settings.plotFpTrace
                obj.fpa.plotTrace();
            end
            
            if obj.settings.plotFpPower
                obj.fpa.plotPowerSpectrum();
            end
            
            if obj.settings.plotFpStats
                obj.fpa.plotStatistics();
            end
            
            if obj.settings.plotFpTrigger
                obj.fpa.plotTrigger();
            end
            
            if obj.settings.plotFpTriggerAverage
                obj.fpa.plotTriggerAverage();
            end
            
            if obj.settings.plotFpAUC
                obj.fpa.plotAUC();
            end
            
            percentile = 0.99;
            grow = 0.50;
            cmap = lines();
            color465 = [0.0000, 0.4470, 0.7410];
            color405 = [0.8500, 0.3250, 0.0980];
            colorDff = color465;
            colorEmg = [1.0000, 0.2500, 0.2500];
            colorEnvelope = [0.0000, 0.0000, 0.0000];
            xlims = obj.fpa.time([1, end]);
            
            if emgRequired && obj.settings.plotEmgTrace
                figureName = 'Raster plots';
                figure('name', figureName);
                
                subplot(3, 1, 1);
                hold('all');
                ylims = limits([obj.fpa.signal; obj.fpa.reference], percentile, grow);
                plotEpochs(epochs, xlims, ylims, cmap, true);
                plot(obj.fpa.time, obj.fpa.signal, 'Color', color465, 'DisplayName', '465');
                plot(obj.fpa.time, obj.fpa.reference, 'Color', color405, 'DisplayName', '405');
                legend('show');
                ylim(ylims);
                
                subplot(3, 1, 2);
                hold('all');
                ylims = limits(obj.fpa.dff, percentile, grow);
                plotEpochs(epochs, xlims, ylims, cmap, false);
                plot(obj.fpa.time, obj.fpa.dff, 'Color', colorDff, 'DisplayName', 'df/f');
                legend('show');
                ylim(ylims);
                
                subplot(3, 1, 3);
                hold('all');
                ylims = limits([emgSignal; emgHigh], percentile, grow);
                plotEpochs(epochs, xlims, ylims, cmap, false);
                plot(emgTime, abs(emgSignal), 'Color', colorEmg, 'DisplayName', 'EMG rectified');
                plot(emgTime, emgHigh, 'Color', colorEnvelope, 'DisplayName', 'EMG envelope');
                legend('show');
                ylim(ylims);
                
                linkaxes(findobj(gcf(), 'type', 'axes'), 'x');
                xlim(xlims);
            end

            if emgRequired && obj.settings.plotXcorr
                figureName = 'Cross-correlation FP to EMG for different behaviors';
                xcLags = round(obj.settings.xcorrLag * obj.settings.resamplingFrequency);
                xcTics = (-xcLags:xcLags) / obj.settings.resamplingFrequency;
                nEpochs = numel(epochs) / 2;
                figure('name', figureName);
                hold('all');
                for e = 1:nEpochs
                    epochName = epochs{2 * e - 1};
                    ranges = epochs{2 * e};
                    mask = time2id(obj.fpa.time, ranges);
                    if isempty(mask)
                        xc = NaN * xcTics;
                    else
                        xc = xcorr(obj.fpa.dff(mask), emgHigh(mask), xcLags, 'normalized');
                    end
                    plot(xcTics, xc, 'DisplayName', epochName);
                end
                title(figureName);
                legend('show');
                xlabel('Lag (s)');
            end
            
            set(obj.h.processButton, 'Enable', true, 'Text', 'Process');
        end
        
        function onGetFile(obj, target)
            switch target
                case obj.h.doricFilenameEdit
                    extensions = {'*.csv'};
                    prompt = 'Select fiber photometry data (Doric)';
                    callback = @obj.onDoricFilenameEdit;
                case obj.h.behaviorFilenameEdit
                    extensions = {'*.csv;*.xlsx;*.tsv', 'BinaryStates (*.csv) or CleverSys (*.xlsx) or BORIS (*.tsv)'};
                    prompt = 'Select behavior data';
                    callback = @obj.onCleverSysFilenameEdit;
                case obj.h.labChartFilenameEdit
                    extensions = {'*.adicht'};
                    prompt = 'Select EMG data (LabChart)';
                    callback = @obj.onLabChartFilenameEdit;
            end
            [file, folder, index] = uigetfile(extensions, prompt, obj.settings.lastFolder);
            if index > 0
                folder = escape(folder);
                obj.settings.lastFolder = folder;
                target.Value = sprintf('%s/%s', folder, file);
                callback();
            end
        end
        
        function onFPBaselineTypeDrop(obj, fpBaselineType)
            obj.h.fpBaselineTypeDrop.Value = fpBaselineType;
            switch fpBaselineType
                case 'exponential decay'
                    %obj.h.fpBaselineLowpassFrequencyEdit.Enable = true;
                    %obj.h.fpBaselineEpochsEdit.Enable = true;
                case 'airPLS'
                    %obj.h.fpBaselineLowpassFrequencyEdit.Enable = false;
                    %obj.h.fpBaselineEpochsEdit.Enable = false;
            end
            obj.saveSettings('fpBaselineType', fpBaselineType);
        end
        
        function onEpochsTypeDrop(obj, epochType)
            obj.h.epochsTypeDrop.Value = epochType;
            obj.h.conditionEpochsPanel.Visible = false;
            obj.h.behaviorPanel.Visible = false;
            switch epochType
                case 'File'
                    obj.h.behaviorPanel.Visible = true;
                case 'Manual'
                    obj.h.conditionEpochsPanel.Visible = true;
            end
            obj.saveSettings('epochsType', epochType);
        end
        
        function onEmgEnvelopeTypeDrop(obj, envelopeType)
            obj.h.emgEnvelopeTypeDrop.Value = envelopeType;
            obj.h.emgEnvelopeLowpassFrequencyEdit.Enable = false;
            obj.h.emgEnvelopeSizeEdit.Enable = false;
            switch envelopeType
                case 'Low-pass'
                    obj.h.emgEnvelopeLowpassFrequencyEdit.Enable = true;
                case 'RMS'
                    obj.h.emgEnvelopeSizeEdit.Enable = true;
            end
            obj.saveSettings('emgEnvelopeType', envelopeType);
        end
        
        function onCleverSysSheetDrop(obj, varargin)
            if numel(varargin) == 1
                sheet = varargin{1};
                obj.h.cleverSysSheetDrop.Value = sheet;
                epochs = loadBehavior(obj.h.behaviorFilenameEdit.Value, sheet);
                entries = epochs(1:2:end);
                obj.behaviorEpochsEntries = entries;
                updateList(obj.h.behaviorEpochsList, obj.behaviorEpochsEntries, obj.settings.behaviorEpochsChoices);
                obj.saveSettings('cleverSysSheet', sheet);
            else
                updateList(obj.h.behaviorEpochsList);
            end
        end
        
        function onCleverSysFilenameEdit(obj)
            target = obj.h.behaviorFilenameEdit;
            filename = target.Value;
            if hasExtension(filename, 'tsv') || hasExtension(filename, 'csv')
                % It's Boris or BinaryStates.
                if exist(obj.h.behaviorFilenameEdit.Value, 'file') == 2
                    epochs = loadBehavior(obj.h.behaviorFilenameEdit.Value);
                    entries = epochs(1:2:end);
                    obj.h.cleverSysSheetDrop.Items = {};
                    obj.onCleverSysSheetDrop();
                    obj.behaviorEpochsEntries = entries;
                    updateList(obj.h.behaviorEpochsList, obj.behaviorEpochsEntries, obj.settings.behaviorEpochsChoices);
                    obj.saveSettings('behaviorFilename', filename);
                    success = true;
                else
                    success = false;
                end
            elseif isempty(filename)
                % User canceled choice.
                obj.h.cleverSysSheetDrop.Items = {};
                obj.onCleverSysSheetDrop();
                obj.saveSettings('behaviorFilename', filename);
                success = true;
            else
                % It's CleverSys.
                try
                    [~, sheets, ~] = xlsfinfo(filename);
                    success = true;
                catch
                    success = false;
                end
                if success
                    obj.h.cleverSysSheetDrop.Items = sheets;
                    if ismember(obj.settings.cleverSysSheet, sheets)
                        sheet = obj.settings.cleverSysSheet;
                    else
                        sheet = sheets{1};
                    end
                    obj.onCleverSysSheetDrop(sheet);
                    obj.saveSettings('behaviorFilename', filename, 'cleverSysSheet', sheet);
                else
                    obj.h.cleverSysSheetDrop.Items = {};
                    obj.onCleverSysSheetDrop();
                end
            end
            setSuccessColor(target, success);
        end
        
        function onLabChartFilenameEdit(obj)
            target = obj.h.labChartFilenameEdit;
            filename = target.Value;
            if isempty(filename)
                updateList(obj.h.labChartChannelsList);
                updateList(obj.h.labChartBlocksList);
                obj.h.labChartCommentsHtml.HTMLSource = '<div></div>';
                obj.saveSettings('labChartFilename', filename);
                success = true;
            else
                try
                    [information, ~, comments] = Aditch.getInformation(filename);
                    success = true;
                catch
                    success = false;
                end
                if success
                    html = '<div style="font-family:monospace">';
                    nComments = numel(comments);
                    for i = 1:2:nComments
                        time = sprintf('%.2f', comments{i});
                        text = comments{i + 1};
                        item = sprintf('<span>[%08s] %s</span>', time, text);
                        html = sprintf('%s%s<br>', html, item);
                    end
                    html = sprintf('%s</div>', html);

                    nNames = numel(information.names);
                    channels = cell(1, nNames);
                    for i = 1:nNames
                        channels{i} = sprintf('#%i "%s"', i, information.names{i});
                    end

                    nBlocks = size(information.period, 2);
                    blocks = cell(1, nBlocks);
                    duration = max(information.samples .* information.period, [], 1);
                    for i = 1:nBlocks
                        from = sprintf('%.2f', information.offset(i));
                        to = sprintf('%.2f', information.offset(i) + duration(i));
                        blocks{i} = sprintf('#%i [%09s .. %09s]', i, from, to);
                    end

                    obj.h.labChartCommentsHtml.HTMLSource = html;
                    obj.labChartChannelsEntries = channels;
                    obj.labChartBlocksEntries = blocks;
                    updateList(obj.h.labChartChannelsList, obj.labChartChannelsEntries, obj.settings.labChartChannelsChoices);
                    updateList(obj.h.labChartBlocksList, obj.labChartBlocksEntries, obj.settings.labChartBlocksChoices);
                    obj.saveSettings('labChartFilename', filename);
                else
                    obj.h.labChartCommentsHtml.HTMLSource = '<div></div>';
                    updateList(obj.h.labChartChannelsList);
                    updateList(obj.h.labChartBlocksList);
                end
            end
            setSuccessColor(target, success);
        end
        
        function onDoricFilenameEdit(obj)
            target = obj.h.doricFilenameEdit;
            filename = target.Value;
            if isempty(filename)
                updateList(obj.h.doricChannelsList);
                obj.saveSettings('doricFilename', filename);
                success = true;
            else
                try
                    [~, channels] = loadDoric(filename);
                    nChannels = numel(channels);
                    for i = 1:nChannels
                        channels{i} = sprintf('#%i %s', i, channels{i});
                    end
                    success = true;
                catch
                    success = false;
                end
                if success
                    obj.doricChannelsEntries = channels;
                    updateList(obj.h.doricChannelsList, obj.doricChannelsEntries, obj.settings.doricChannelsChoices);
                    obj.saveSettings('doricFilename', filename);
                else
                    updateList(obj.h.doricChannelsList);
                end
            end
            setSuccessColor(target, success);
        end
        
        function onBaselineEpochsEdit(obj)
            target = obj.h.fpBaselineEpochsEdit;
            text = target.Value;
            % Turn contents into a matrix.
            text = ['[', regexprep(text, '^\s*\[|]\s*$', ''), ']'];
            try
                epochRanges = eval(text);
                target.Value = text;
                success = validateEpochRanges(epochRanges);
            catch
                success = false;
            end
            setSuccessColor(target, success);
            if success
                obj.saveSettings('fpBaselineEpochs', epochRanges, 'fpBaselineEpochsText', text);
            end
        end
        
        function onArtifactEpochsEdit(obj)
            target = obj.h.fpArtifactEpochsEdit;
            text = target.Value;
            % Turn contents into a matrix.
            text = ['[', regexprep(text, '^\s*\[|]\s*$', ''), ']'];
            try
                epochRanges = eval(text);
                target.Value = text;
                success = validateEpochRanges(epochRanges);
            catch
                success = false;
            end
            setSuccessColor(target, success);
            if success
                obj.saveSettings('fpArtifactEpochs', epochRanges, 'fpArtifactEpochsText', text);
            end
        end
        
        function onConditionEpochsEdit(obj)
            target = obj.h.conditionEpochsEdit;
            % Merge lines.
            text = target.Value;
            text = strjoin(text, ' ');
            % Turn contents into a cell.
            text = ['{', regexprep(text, '^\s*\{|}\s*$', ''), '}'];
            try
                epochs = eval(text);
                target.Value = text;
                success = validateEpochs(epochs);
            catch
                success = false;
            end
            setSuccessColor(target, success);
            if success
                obj.saveSettings('conditionEpochs', epochs, 'conditionEpochsText', text);
            end
        end
        
        function setChoice(obj, list, choice)
            namesShowing = list.Items;
            index = ismember(namesShowing, list.Value);
            switch list
                case obj.h.doricChannelsList
                    obj.settings.doricChannelsChoices = updateChoices(obj.settings.doricChannelsChoices, obj.doricChannelsEntries(index), choice);
                    updateList(list, obj.doricChannelsEntries, obj.settings.doricChannelsChoices);
                    obj.saveSettings('doricChannelsChoices', obj.settings.doricChannelsChoices);
                case obj.h.behaviorEpochsList
                    obj.settings.behaviorEpochsChoices = updateChoices(obj.settings.behaviorEpochsChoices, obj.behaviorEpochsEntries(index), choice);
                    updateList(list, obj.behaviorEpochsEntries, obj.settings.behaviorEpochsChoices);
                    obj.saveSettings('behaviorEpochsChoices', obj.settings.behaviorEpochsChoices);
                case obj.h.labChartChannelsList
                    obj.settings.labChartChannelsChoices = updateChoices(obj.settings.labChartChannelsChoices, obj.labChartChannelsEntries(index), choice);
                    updateList(list, obj.labChartChannelsEntries, obj.settings.labChartChannelsChoices);
                    obj.saveSettings('labChartChannelsChoices', obj.settings.labChartChannelsChoices);
                case obj.h.labChartBlocksList
                    obj.settings.labChartBlocksChoices = updateChoices(obj.settings.labChartBlocksChoices, obj.labChartBlocksEntries(index), choice);
                    updateList(list, obj.labChartBlocksEntries, obj.settings.labChartBlocksChoices);
                    obj.saveSettings('labChartBlocksChoices', obj.settings.labChartBlocksChoices);
            end
        end
        
        function contextMenu = addMenu(obj, target, options, callback)
            contextMenu = uicontextmenu(obj.h.control);
            target.ContextMenu = contextMenu;
            display = options(1:2:end);
            aliases = options(2:2:end);
            nOptions = numel(display);
            for i = 1:nOptions
                menu = uimenu(contextMenu);
                menu.MenuSelectedFcn = @(~, ~)callback(target, aliases{i});
                menu.Text = display{i};
            end
        end
        
        function onClose(~)
            delete(obj.h.control);
        end
    end
    
    methods
        function export(obj, varargin)
            if numel(varargin) == 1
                prefix = varargin;
            else
                prefix = fullfile(obj.settings.lastFolder, 'output');
            end
            obj.fpa.export(prefix);
        end
    end
    
    methods (Static)
        function options = enumerate(type, varargin)
            switch type
                case 'epochsType'
                    options = {'File', 'Manual'};
                case 'emgBandpassType'
                    options = {'Zero-phase', 'Butter', 'Window'};
                case 'emgEnvelopeType'
                    options = {'Low-pass', 'RMS'};
                case 'fpBaselineType'
                    options = {'exponential decay', 'airPLS'};
                case 'fpFitReference'
                    options = {'true', 'false'};
                case 'fpThresholdingFunctions'
                    options = {'median:mad', 'mean:std'};
                case 'fpNormalizationType'
                    options = {'z-score', 'df/f', 'median:mad', 'median:std', 'mean:mad'};
                otherwise
                    options = {};
            end
        end
    end
end

function settings = getDefaults()
    % Internal.
    settings.lastFolder = fullfile(getHome(), 'Documents');
    
    % General.
    settings.epochsType = 'Manual';
            
    % Frequency to resample EMG and FP.
    settings.resamplingFrequency = 100;
    % Maximum cross-correlation lag between EMG and FP.
    settings.xcorrLag = 5;
    
    % EMG.
    % How to filter EMG.
    settings.emgBandpassType = 'Zero-phase';
    % How to calculate envelope.
    settings.emgEnvelopeType = 'RMS';
    % Frequency of the envelope filter when choice is low-pass.
    settings.emgEnvelopeLowpassFrequency = 1;
    % Duration of the envelope window when choice is rms.
    settings.emgEnvelopeSize = 1;
    % Burden et al 2003. Muscles: vastus lateralis, vastus medialis, semitendinosus, biceps femoris.
    %   Under 20Hz: noise.
    %   Under 450Hz: still clearly identifiable rapid on-off bursts of the EMG.
    settings.emgBandpassFrequencyLow = 20;
    settings.emgBandpassFrequencyHigh = 120;
    
    % FP.
    settings.fpFilename = '';
    settings.fpBaselineType = 'exponential decay';
    settings.fpBaselineEpochs = [];
    settings.fpBaselineEpochsText = '';
    settings.fpBaselineLowpassFrequency = 0.1;
    settings.fpFitReference = 'true';
    settings.fpArtifactEpochs = [];
    settings.fpArtifactEpochsText = '';
    settings.fpLowpassFrequency = 2;
    settings.fpPeaksLowpassFrequency = 0.5;
    settings.fpThresholdingFunctions = 'median:mad';
    settings.fpThresholdingFactor = 2.91;
    settings.fpNormalizationType = 'median:mad';
    
    % Epochs.
    settings.conditionEpochs = {};
    settings.conditionEpochsText = '';
    
    % CleverSys epochs.
    settings.behaviorFilename = '';
    settings.cleverSysSheet = '';
    settings.behaviorEpochsChoices = cell(0, 2);
    
    % Doric.
    settings.doricFilename = '';
    settings.doricChannelsChoices = cell(0, 2);
    
    % LabChart.
    settings.labChartFilename = '';
    settings.labChartChannelsChoices = cell(0, 2);
    settings.labChartBlocksChoices = cell(0, 2);
    
    % Plots.
    settings.plotFpTrace = true;
    settings.plotFpPower = true;
    settings.plotFpStats = true;
    settings.plotFpTrigger = true;
    settings.plotFpTriggerAverage = true;
    settings.plotFpAUC = true;
    settings.plotEmgTrace = true;
    settings.plotXcorr = true;
end

function set(objs, varargin)
    nArgs = numel(varargin);
    nObjs = numel(objs);
    for i = 1:nObjs
        for j = 1:2:nArgs
            key = varargin{j};
            value = varargin{j + 1};
            objs(i).(key) = value;
        end
    end
end
        
function updateList(list, varargin)
    if numel(varargin) == 0
        [entries, choices] = deal({}, cell(0, 2));
    else
        [entries, choices] = deal(varargin{:});
    end
    labels = cell(size(entries));
    charCounts = cellfun(@numel, entries);
    tabCounts = max(charCounts) - charCounts + 2;
    nEntries = numel(entries);
    keys = choices(:, 1);
    values = choices(:, 2);
    for i = 1:nEntries
        entry = entries{i};
        k = ismember(keys, entry);
        if any(k)
            value = values{k};
        else
            value = 'Unset';
        end
        if isequal(lower(value), 'unset')
            label = entry;
        else
            label = sprintf('%s%*s==> %s', entry, tabCounts(i), ' ', value);
        end
        labels{i} = label;
    end
    index = ismember(list.Items, list.Value);
    index(nEntries + 1:end) = [];
    list.Items = labels;
    list.Value = labels(index);
end

function choicePairs = updateChoices(choicePairs, updateKeys, updateValue)
    choiceKeys = choicePairs(:, 1);
    choiceValues = choicePairs(:, 2);
    existingKeyMask = ismember(choiceKeys, updateKeys);
    [~, addKeyIds] = setdiff(updateKeys, choiceKeys);
    
    [choiceValues{existingKeyMask}] = deal(updateValue);
    addKeys = updateKeys(addKeyIds);
    nAdd = numel(addKeys);
    choiceKeys(end + 1:end + nAdd) = addKeys;
    [choiceValues{end + 1:end + nAdd}] = deal(updateValue);
    choicePairs = [choiceKeys(:), choiceValues(:)];
end

function indexes = getIndex(entries, choicePairs, value)
    choiceKeys = choicePairs(:, 1);
    choiceValues = choicePairs(:, 2);
    entryKeys = entries;
    entryValues = repmat({''}, size(entryKeys));
    [~, e, c] = intersect(entryKeys, choiceKeys);
    entryValues(e) = choiceValues(c);
    indexes = find(ismember(entryValues, value));
end

function setSuccessColor(target, success)
    if success
        target.BackgroundColor = [1.0, 1.0, 1.0];
    else
        target.BackgroundColor = [1.0, 0.5, 0.5];
    end
end

function text = escape(text)
    text = strrep(text, '\', '/');
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

function [time2, varargout] = alignTimestamps(time, dt, varargin)
    % Shift points so that first time point is a multiple of dt.
    timestamp0 = ceil(time(1) / dt) * dt;
    time2 = colon(timestamp0, median(diff(time)), time(end))';
    varargout = cellfun(@(data) interp1(time(:), data(:), time2), varargin, 'UniformOutput', false);
end

function ylims = limits(x, percentile, grow)
    ylims = [prctile(x, 100 * (1 - percentile)), prctile(x, 100 * percentile)];
    delta = diff(ylims) * grow;
    ylims = [ylims(1) - delta, ylims(2) + delta];
end

function errorDialog(messages)
    errordlg(messages, 'Processing error', 'modal');
end

function test = hasExtension(filename, extension)
    [~, ~, ext] = fileparts(filename);
    test = lower(ext) == "." + extension;
end

function varargout = loadBehavior(varargin)
    filename = varargin{1};
    if hasExtension(filename, 'tsv')
        [varargout{1:nargout}] = loadBoris(filename);
    elseif hasExtension(filename, 'csv')
        [varargout{1:nargout}] = loadBinaryStates(filename);
    else
        [varargout{1:nargout}] = loadCleverSys(varargin{:});
    end
end

function folder = getHome()
    if ispc
        folder = getenv('USERPROFILE');
    else
        folder = char(java.lang.System.getProperty('user.home'));
    end
end

function [data, channels] = loadDoric(filename)
    data = readtable(filename);
    channels = data.Properties.VariableNames;
    data = data{:, :};
    % Remove NaN.
    k = any(isnan(data), 2);
    data = data(~k, :);
end

%#ok<*PROPLC>