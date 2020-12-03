% GUI for processing EMG and FP data.
% see FPA.

% 2020-11-01. Leonardo Molina.
% 2020-12-02. Last modified.
classdef GUI < handle
    properties % (Access = private)
        doricChannelsEntries
        cleverSysEpochsEntries
        labChartChannelsEntries
        labChartBlocksEntries
        conditionEpochs
        settingsFilename
        cleverSysIsCSV
        settings
        cache
        h
    end
    
    methods
        function obj = GUI(varargin)
            addpath('../../');
            addpath('../../common');
            
            obj.doricChannelsEntries = {};
            obj.cleverSysEpochsEntries = {};
            obj.labChartChannelsEntries = {};
            obj.labChartBlocksEntries = {};
            obj.conditionEpochs = [-Inf, Inf];
            obj.cache.doricFilename = '<empty>';
            obj.cache.labChartFilename = '<empty>';
            obj.cache.cleverSysFilename = '<empty>';
            
            % Calculate GUI dimensions.
            screenSize = get(0, 'ScreenSize');
            screenSize = screenSize(3:4);
            targetSize = screenSize([1, 2]) .* [0.75, 0.75];
            uiHeight = 27;
            nRows = floor(targetSize(2) / uiHeight) - 2;
            
            % Create main figure.
            obj.h.control = uifigure('Name', 'EMG * FP analysis', 'MenuBar', 'none', 'NumberTitle', 'off', 'ToolBar', 'none', 'WindowState', 'maximized'); % , 'CloseRequestFcn', @(~, ~)obj.onClose
            obj.h.control.Position = [(screenSize(1) - targetSize(1)) / 2, (screenSize(2) - targetSize(2)) / 2, targetSize(1), targetSize(2)];
            
            % 3x3 layout.
            rowCount = floor(nRows / 3);
            mainLayout = uigridlayout(obj.h.control);
            mainLayout.RowHeight = repmat({'1x'}, 1, 3);
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
            label = uilabel(panel, 'Text', 'Behavior duration', 'HorizontalAlignment', 'right');
            label.Layout.Row = 3;
            label.Layout.Column = 1;
            label = uilabel(panel, 'Text', 'Cross-correlation lag', 'HorizontalAlignment', 'right');
            label.Layout.Row = 4;
            label.Layout.Column = 1;
            
            % Epochs type.
            h = uidropdown(panel, 'Items', {'CleverSys', 'Manual'}, 'ValueChangedFcn', @(h, ~)obj.onEpochsTypeDrop(h.Value));
            h.Layout.Row = 1;
            h.Layout.Column = 2;
            obj.h.epochsTypeDrop = h;
            h = uieditfield(panel, 'numeric', 'Limits', [10, Inf], 'ValueChangedFcn', @(h, ~)obj.saveSettings('resamplingFrequency', h.Value));
            h.Layout.Row = 2;
            h.Layout.Column = 2;
            obj.h.resamplingFrequencyEdit = h;
            h = uieditfield(panel, 'numeric', 'Limits', [0, Inf], 'ValueChangedFcn', @(h, ~)obj.saveSettings('behaviorDuration', h.Value));
            h.Layout.Row = 3;
            h.Layout.Column = 2;
            obj.h.behaviorDurationEdit = h;
            h = uieditfield(panel, 'numeric', 'Limits', [1, Inf], 'ValueChangedFcn', @(h, ~)obj.saveSettings('xcorrLag', h.Value));
            h.Layout.Row = 4;
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
            
            h = uidropdown(panel, 'Items', {'Zero-phase', 'Butter', 'Window'}, 'ValueChangedFcn', @(h, ~)obj.saveSettings('emgBandpassType', h.Value));
            h.Layout.Row = 1;
            h.Layout.Column = [2, 3];
            obj.h.emgBandpassTypeDrop = h;
            h = uidropdown(panel, 'Items', {'Low-pass', 'RMS'}, 'ValueChangedFcn', @(h, ~)obj.onEmgEnvelopeTypeDrop(h.Value));
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
            panel.ColumnWidth = {'4x', '1x', '1x'};
            
            % FP settings.
            label.Layout.Column = 1;
            label = uilabel(panel, 'Text', 'Bleaching epochs', 'HorizontalAlignment', 'right');
            label.Layout.Row = 1;
            label.Layout.Column = 1;
            label = uilabel(panel, 'Text', 'Peaks bandpass frequency', 'HorizontalAlignment', 'right');
            label.Layout.Row = 2;
            label.Layout.Column = 1;
            label = uilabel(panel, 'Text', 'df/f lowpass frequency', 'HorizontalAlignment', 'right');
            label.Layout.Row = 3;
            label.Layout.Column = 1;
            label = uilabel(panel, 'Text', 'Bleaching lowpass frequency', 'HorizontalAlignment', 'right');
            label.Layout.Row = 4;
            label.Layout.Column = 1;
            label = uilabel(panel, 'Text', 'Thresholding function', 'HorizontalAlignment', 'right');
            label.Layout.Row = 5;
            label.Layout.Column = 1;
            label = uilabel(panel, 'Text', 'Thresholding factor', 'HorizontalAlignment', 'right');
            label.Layout.Row = 6;
            label.Layout.Column = 1;
            label = uilabel(panel, 'Text', 'Normalization type', 'HorizontalAlignment', 'right');
            label.Layout.Row = 7;
            label.Layout.Column = 1;
            label = uilabel(panel, 'Text', 'Window size', 'HorizontalAlignment', 'right');
            label.Layout.Row = 8;
            label.Layout.Column = 1;
            
            h = uieditfield(panel, 'ValueChangedFcn', @(~, ~)obj.onBleachingEpochsEdit);
            h.Layout.Row = 1;
            h.Layout.Column = [2, 3];
            obj.h.fpBleachingEpochsEdit = h;
            h = uieditfield(panel, 'numeric', 'Limits', [1e-3, Inf], 'ValueChangedFcn', @(h, ~)obj.saveSettings('fpPeaksBandpassFrequencyLow', h.Value));
            h.Layout.Row = 2;
            h.Layout.Column = 2;
            obj.h.fpPeaksBandpassFrequencyLowEdit = h;
            h = uieditfield(panel, 'numeric', 'Limits', [1e-3, Inf], 'ValueChangedFcn', @(h, ~)obj.saveSettings('fpPeaksBandpassFrequencyHigh', h.Value));
            h.Layout.Row = 2;
            h.Layout.Column = 3;
            obj.h.fpPeaksBandpassFrequencyHighEdit = h;
            h = uieditfield(panel, 'numeric', 'Limits', [1e-3, Inf], 'ValueChangedFcn', @(h, ~)obj.saveSettings('fpDffLowpassFrequency', h.Value));
            h.Layout.Row = 3;
            h.Layout.Column = [2, 3];
            obj.h.fpDffLowpassFrequencyEdit = h;
            h = uieditfield(panel, 'numeric', 'Limits', [1e-3, Inf], 'ValueChangedFcn', @(h, ~)obj.saveSettings('fpBleachingLowpassFrequency', h.Value));
            h.Layout.Row = 4;
            h.Layout.Column = [2, 3];
            obj.h.fpBleachingLowpassFrequencyEdit = h;
            h = uidropdown(panel, 'Items', {'mad', 'std'}, 'ValueChangedFcn', @(h, ~)obj.saveSettings('fpThresholdingFunction', h.Value));
            h.Layout.Row = 5;
            h.Layout.Column = [2, 3];
            obj.h.fpThresholdingFunctionDrop = h;
            h = uieditfield(panel, 'numeric', 'Limits', [1e-3, Inf], 'ValueChangedFcn', @(h, ~)obj.saveSettings('fpThresholdingFactor', h.Value));
            h.Layout.Row = 6;
            h.Layout.Column = [2, 3];
            obj.h.fpThresholdingFactorEdit = h;
            h = uidropdown(panel, 'Items', {'z-score', 'df/f'}, 'ValueChangedFcn', @(h, ~)obj.saveSettings('fpNormalizationType', h.Value));
            h.Layout.Row = 7;
            h.Layout.Column = [2, 3];
            obj.h.fpNormalizationTypeDrop = h;
            h = uieditfield(panel, 'numeric', 'Limits', [1e-3, Inf], 'ValueChangedFcn', @(h, ~)obj.saveSettings('fpWindowSize', h.Value));
            h.Layout.Row = 8;
            h.Layout.Column = [2, 3];
            obj.h.fpWindowSizeEdit = h;
            
            % Epochs panel.
            panel = uipanel(mainLayout, 'Title', 'Manual Epoch definitions');
            obj.h.conditionEpochsPanel = panel;
            panel = uigridlayout(panel, [1, 1]);
            panel.RowHeight = {'1x'};
            panel.ColumnWidth = {'1x'};
            h = uitextarea(panel, 'ValueChangedFcn', @(~, ~)obj.onConditionEpochsEdit);
            obj.h.conditionEpochsEdit = h;
            
            % CleverSys epochs panel.
            panel = uipanel(mainLayout, 'Title', 'CleverSys Epoch definitions');
            panel.Layout = obj.h.conditionEpochsPanel.Layout;
            obj.h.cleverSysPanel = panel;
            panel = uigridlayout(panel, [1, 1]);
            panel.RowHeight = repmat({'1x'}, 1, rowCount);
            panel.ColumnWidth = {'5x', '1x'};
            
            % CleverSys filename.
            h = uieditfield(panel, 'ValueChangedFcn', @(~, ~)obj.onCleverSysFilenameEdit);
            h.Layout.Row = 1;
            h.Layout.Column = 1;
            obj.h.cleverSysFilenameEdit = h;
            h = uibutton(panel, 'Text', 'Select', 'ButtonPushedFcn', @(~, ~)obj.onGetFile(obj.h.cleverSysFilenameEdit));
            h.Layout.Row = 1;
            h.Layout.Column = 2;
            obj.h.cleverSysButton = h;
            
            % CleverSys sheet dropdown and list.
            h = uidropdown(panel, 'Items', {}, 'ValueChangedFcn', @(h, ~)obj.onCleverSysSheetDrop(h.Value));
            h.Layout.Row = 2;
            h.Layout.Column = [1, 2];
            obj.h.cleverSysSheetDrop = h;
            h = uilistbox(panel, 'Items', {}, 'Multiselect', true, 'FontName', 'Monospaced');
            h.Layout.Row = [3, rowCount];
            h.Layout.Column = [1, 2];
            obj.addMenu(h, {'Use', 'Use', 'Unset', 'Unset'}, @obj.setChoice);
            obj.h.cleverSysEpochsList = h;
            
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
            obj.addMenu(h, {'Set as signal', 'Signal', 'Set as reference', 'Reference', 'Set as camera', 'Camera', 'Unset', 'Unset'}, @obj.setChoice);
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
            panel.RowHeight = repmat({'1x'}, 1, rowCount);
            panel.ColumnWidth = {'1x'};
            
            savePlotSelection = @()obj.saveSettings('plotFpTrace', obj.h.fpPlots(1).Value, 'plotFpPower', obj.h.fpPlots(2).Value, 'plotFpStats', obj.h.fpPlots(3).Value, 'plotFpTrigger', obj.h.fpPlots(4).Value, 'plotEmgTrace', obj.h.emgPlots(1).Value, 'plotXcorr', obj.h.generalPlots(1).Value);
            obj.h.fpPlots(1) = uicheckbox(panel, 'Text', 'df/f trace and peaks', 'ValueChangedFcn', @(~, ~)savePlotSelection());
            obj.h.fpPlots(2) = uicheckbox(panel, 'Text', 'df/f power spectrum', 'ValueChangedFcn', @(~, ~)savePlotSelection());
            obj.h.fpPlots(3) = uicheckbox(panel, 'Text', 'df/f stats boxplot', 'ValueChangedFcn', @(~, ~)savePlotSelection());
            obj.h.fpPlots(4) = uicheckbox(panel, 'Text', 'df/f triggered average', 'ValueChangedFcn', @(~, ~)savePlotSelection());
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
                obj.settingsFilename = fullfile(getenv('USERPROFILE'), 'Documents', 'MATLAB', 'csmopto-settings.mat');
                if exist(obj.settingsFilename, 'file') == 2
                    settings = load(obj.settingsFilename);
                    settings = settings.settings;
                    fprintf('Loaded settings from "%s"\n', obj.settingsFilename);
                else
                    settings = defaults();
                end
            end
            obj.applySettings(settings);
        end
    end
    
    methods
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
        
        function applySettings(obj, settings)
            obj.settings = settings;
            
            obj.h.epochsTypeDrop.Value = settings.epochsType;
            obj.h.resamplingFrequencyEdit.Value = settings.resamplingFrequency;
            obj.h.behaviorDurationEdit.Value = settings.behaviorDuration;
            obj.h.xcorrLagEdit.Value = settings.xcorrLag;
            
            obj.h.emgBandpassTypeDrop.Value = settings.emgBandpassType;
            obj.h.emgEnvelopeTypeDrop.Value = settings.emgEnvelopeType;
            obj.h.emgEnvelopeLowpassFrequencyEdit.Value = settings.emgEnvelopeLowpassFrequency;
            obj.h.emgEnvelopeSizeEdit.Value = settings.emgEnvelopeSize;
            obj.h.emgBandpassFrequencyLowEdit.Value = settings.emgBandpassFrequencyLow;
            obj.h.emgBandpassFrequencyHighEdit.Value = settings.emgBandpassFrequencyHigh;
            
            obj.h.fpBleachingEpochsEdit.Value = settings.fpBleachingEpochs;
            obj.h.fpPeaksBandpassFrequencyLowEdit.Value = settings.fpPeaksBandpassFrequencyLow;
            obj.h.fpPeaksBandpassFrequencyHighEdit.Value = settings.fpPeaksBandpassFrequencyHigh;
            obj.h.fpDffLowpassFrequencyEdit.Value = settings.fpDffLowpassFrequency;
            obj.h.fpBleachingLowpassFrequencyEdit.Value = settings.fpBleachingLowpassFrequency;
            obj.h.fpThresholdingFunctionDrop.Value = settings.fpThresholdingFunction;
            obj.h.fpThresholdingFactorEdit.Value = settings.fpThresholdingFactor;
            obj.h.fpNormalizationTypeDrop.Value = settings.fpNormalizationType;
            obj.h.fpWindowSizeEdit.Value = settings.fpWindowSize;
            
            obj.h.conditionEpochsEdit.Value = settings.conditionEpochs;
            
            obj.h.cleverSysFilenameEdit.Value = settings.cleverSysFilename;
            obj.h.doricFilenameEdit.Value = settings.doricFilename;
            obj.h.labChartFilenameEdit.Value = settings.labChartFilename;
            
            obj.h.fpPlots(1).Value = settings.plotFpTrace;
            obj.h.fpPlots(2).Value = settings.plotFpPower;
            obj.h.fpPlots(3).Value = settings.plotFpStats;
            obj.h.fpPlots(4).Value = settings.plotFpTrigger;
            obj.h.emgPlots(1).Value = settings.plotEmgTrace;
            obj.h.generalPlots(1).Value = settings.plotXcorr;
            
            % Cascade calls.
            obj.onDoricFilenameEdit();
            obj.onCleverSysFilenameEdit();
            obj.onLabChartFilenameEdit();
            obj.onEpochsTypeDrop(settings.epochsType);
            obj.onEmgEnvelopeTypeDrop(settings.emgEnvelopeType);
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
            
            set(obj.h.processButton, 'Enable', false, 'Text', 'Processing...');
            drawnow();
            messages = {};

            % Fiber-photometry settings.
            fp = struct();
            fp.signalChannel = find(cellfun(@(option) isequal(option, 'Signal'), obj.settings.doricChannelsChoices));
            fp.referenceChannel = find(cellfun(@(option) isequal(option, 'Reference'), obj.settings.doricChannelsChoices));
            fp.cameraChannel = find(cellfun(@(option) isequal(option, 'Camera'), obj.settings.doricChannelsChoices));
            if numel(fp.signalChannel) ~= 1
                messages = cat(2, messages, 'Select a single Signal channel in Doric Channels.');
            end
            if numel(fp.referenceChannel) ~= 1
                messages = cat(2, messages, 'Select a single Reference channel in Doric Channels.');
            end
            if numel(fp.cameraChannel) ~= 1
                messages = cat(2, messages, 'Select a single Camera channel in Doric Channels.');
            end
            if numel(messages) > 0
                errorDialog(messages);
                set(obj.h.processButton, 'Enable', true, 'Text', 'Process');
                return;
            end
            
            fp.bleachingEpochs = eval(obj.settings.fpBleachingEpochs);
            fp.peaksBandpassFrequency = [obj.settings.fpPeaksBandpassFrequencyLow, obj.settings.fpPeaksBandpassFrequencyHigh];
            fp.dffLowpassFrequency = obj.settings.fpDffLowpassFrequency;
            fp.bleachingLowpassFrequency = obj.settings.fpBleachingLowpassFrequency;
            fp.resamplingFrequency = obj.settings.resamplingFrequency;
            switch obj.settings.fpThresholdingFunction
                case 'mad'
                    fp.thresholdingFunction = @mad;
                case 'std'
                    fp.thresholdingFunction = @std;
            end
            fp.thresholdFactor = obj.settings.fpThresholdingFactor;
            switch obj.settings.fpNormalizationType
                case 'z-score'
                    fp.f0Function = @movmean;
                    fp.f1Function = @movstd;
                case 'df/f'
                    fp.f0Function = @movmean;
                    fp.f1Function = @movmean;
            end
            fp.f0Window = obj.settings.fpWindowSize;
            fp.f1Window = obj.settings.fpWindowSize;
            
            % LabChart settings.
            labChart = struct();
            labChart.emgChannel = find(cellfun(@(option) isequal(option, 'EMG'), obj.settings.labChartChannelsChoices));
            labChart.block = find(cellfun(@(option) isequal(option, 'Use'), obj.settings.labChartBlocksChoices));
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
            
            fprintf('Loading Doric data ... ');
            if isequal(obj.settings.doricFilename, obj.cache.doricFilename)
                fprintf('reused cache.\n');
                fpData = obj.cache.fpData;
            else
                fprintf('\n');
                fpData = loadData(obj.settings.doricFilename);
                obj.cache.fpData = fpData;
                obj.cache.doricFilename = obj.settings.doricFilename;
            end
            
            fpTime = fpData(:, 1);
            fpSignal = fpData(:, fp.signalChannel);
            fpReference = fpData(:, fp.referenceChannel);
            
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
            
            fprintf('Processing data ...\n');
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
            
            switch obj.settings.epochsType
                case 'Manual'
                    epochs = eval(obj.settings.conditionEpochs);
                case 'CleverSys'
                    fprintf('Loading CleverSys data ... ');
                    if isequal(obj.settings.cleverSysFilename, obj.cache.cleverSysFilename) && isequal(obj.settings.cleverSysSheet, obj.cache.cleverSysSheet)
                        fprintf('reused cache.\n');
                        epochs = obj.cache.cleverSysEpochs;
                    else
                        fprintf('\n');
                        epochs = loadCleverSys(obj.settings.cleverSysFilename, obj.settings.cleverSysSheet);
                        obj.cache.cleverSysEpochs = epochs;
                        obj.cache.cleverSysFilename = obj.settings.cleverSysFilename;
                        obj.cache.cleverSysSheet = obj.settings.cleverSysSheet;
                    end
                    k = find(cellfun(@(option) isequal(option, 'Use'), obj.settings.cleverSysEpochsChoices));
                    k = sort([2 * k, 2 * k - 1]);
                    epochs = epochs(k);
                    
                    fpCameraData = fpData(:, fp.cameraChannel);
                    behaviorStart = fpTime(find(fpCameraData, 1));
                    epochs(2:2:end) = cellfun(@(epoch) epoch + behaviorStart, epochs(2:2:end), 'UniformOutput', false);
                    behaviorEnd = max(cat(2, epochs{2:2:end}, behaviorStart + obj.settings.behaviorDuration));
                    preBaseline = [-Inf, behaviorStart];
                    postBaseline = [behaviorEnd, Inf];
                    epochs = ['pre-baseline', preBaseline, epochs, 'post-baseline', postBaseline];
            end
            
            fp.conditionEpochs = epochs;
            fp.plot = {'trace', 'power', 'stats', 'trigger'};
            fp.plot = fp.plot([obj.settings.plotFpTrace, obj.settings.plotFpPower, obj.settings.plotFpStats, obj.settings.plotFpTrigger]);
            [fpTime, fpSignal, fpReference] = alignTimestamps(fpTime, 1 / obj.settings.resamplingFrequency, fpSignal, fpReference);
            fpa = FPA(fpTime, fpSignal, fpReference, fp);
            cellfun(@warning, fpa.warnings);
            
            [emgTime, emgSignal] = alignTimestamps(emgTime, 1 / obj.settings.resamplingFrequency, emgSignal);
            sourceFrequency = 1 / median(diff(emgTime));
            
            % Bandpass emg to remove artifacts.
            switch obj.settings.emgBandpassType
                case 'zero-phase'
                    % High-order, butter filter, without phase shift.
                    bandpassFilter = designfilt('bandpassiir', 'HalfPowerFrequency1', obj.settings.emgBandpassFrequency(1), 'HalfPowerFrequency2', obj.settings.emgBandpassFrequency(2), 'SampleRate', sourceFrequency, 'DesignMethod', 'butter', 'FilterOrder', 12);
                    emgSignal = filtfilt(bandpassFilter, emgSignal);
                case 'butter'
                    % Standard butter filter.
                    % http://www1.udel.edu/biology/rosewc/kaap686/notes/EMG%20analysis.pdf
                    fn = sourceFrequency / 2;
                    fl = obj.settings.emgBandpassFrequency(1) / fn;
                    fh = obj.settings.emgBandpassFrequency(2) / fn;
                    [b, a] = butter(4, [fl, fh], 'bandpass');
                    emgSignal = filter(b, a, emgSignal);
                case 'window'
                    % MATLAB's default filter.
                    bandpassFilter = designfilt('bandpassfir', 'CutoffFrequency1', obj.settings.emgBandpassFrequency(1), 'CutoffFrequency2', obj.settings.emgBandpassFrequency(2), 'SampleRate', sourceFrequency, 'DesignMethod', 'window', 'FilterOrder', 12);
                    emgSignal = filtfilt(bandpassFilter, emgSignal);
            end

            % Notch filter to remove 60Hz noise.
            notchFilter = designfilt('bandstopiir', 'FilterOrder', 2, 'HalfPowerFrequency1', 59, 'HalfPowerFrequency2', 61, 'DesignMethod', 'butter', 'SampleRate', sourceFrequency);
            emgSignal = filtfilt(notchFilter, emgSignal);
            % Detrend.
            emgSignal = detrend(emgSignal);
            % Make data have the same sampling rate.
            if obj.settings.resamplingFrequency < sourceFrequency
                [p, q] = rat(obj.settings.resamplingFrequency / sourceFrequency);
                [emgSignal, emgTime] = resample(emgSignal, emgTime, obj.settings.resamplingFrequency, p, q);
            else
                warning('[emg] Cannot resample to frequencies higher than the source frequency (%.2f Hz).', sourceFrequency);
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
            
            percentile = 0.99;
            grow = 0.50;
            cmap = lines();
            color465 = [0.0000, 0.4470, 0.7410];
            color405 = [0.8500, 0.3250, 0.0980];
            colorDff = color465;
            colorEmg = [1.0000, 0.2500, 0.2500];
            colorEnvelope = [0.0000, 0.0000, 0.0000];
            xlims = fpa.time([1, end]);
            
            if obj.settings.plotEmgTrace
                figureName = 'Raster plots';
                figure('name', figureName);
                
                subplot(3, 1, 1);
                hold('all');
                ylims = limits([fpa.signal; fpa.reference], percentile, grow);
                plotEpochs(epochs, xlims, ylims, cmap, true);
                plot(fpa.time, fpa.signal, 'Color', color465, 'DisplayName', '465');
                plot(fpa.time, fpa.reference, 'Color', color405, 'DisplayName', '405');
                legend('show');
                ylim(ylims);
                
                subplot(3, 1, 2);
                hold('all');
                ylims = limits(fpa.dff, percentile, grow);
                plotEpochs(epochs, xlims, ylims, cmap, false);
                plot(fpa.time, fpa.dff, 'Color', colorDff, 'DisplayName', 'df/f');
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

            if obj.settings.plotXcorr
                figureName = 'Cross-correlation FP to EMG for different behaviors';
                xcLags = round(obj.settings.xcorrLag * obj.settings.resamplingFrequency);
                xcTics = (-xcLags:xcLags) / obj.settings.resamplingFrequency;
                nEpochs = numel(epochs) / 2;
                figure('name', figureName);
                hold('all');
                for e = 1:nEpochs
                    epochName = epochs{2 * e - 1};
                    ranges = epochs{2 * e};
                    mask = time2id(fpa.time, ranges);
                    if isempty(mask)
                        xc = NaN * xcTics;
                    else
                        xc = xcorr(fpa.dff(mask), emgHigh(mask), xcLags, 'normalized');
                    end
                    plot(xcTics, xc, 'DisplayName', epochName);
                end
                title(figureName);
                legend('show');
                xlabel('Lag (s)');
            end
            
            set(obj.h.processButton, 'Enable', true, 'Text', 'Process');
        end
        
        function onEpochsTypeDrop(obj, epochType)
            obj.h.conditionEpochsPanel.Visible = false;
            obj.h.cleverSysPanel.Visible = false;
            switch epochType
                case 'CleverSys'
                    obj.h.cleverSysPanel.Visible = true;
                    obj.h.behaviorDurationEdit.Enable = true;
                case 'Manual'
                    obj.h.conditionEpochsPanel.Visible = true;
                    obj.h.behaviorDurationEdit.Enable = false;
            end
            obj.saveSettings('epochsType', epochType);
        end
        
        function onEmgEnvelopeTypeDrop(obj, envelopeType)
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
        
        function onGetFile(obj, target)
            switch target
                case obj.h.doricFilenameEdit
                    extensions = {'*.csv'};
                    prompt = 'Select fiber photometry data (Doric)';
                    callback = @obj.onDoricFilenameEdit;
                case obj.h.cleverSysFilenameEdit
                    extensions = {'*.xlsx'};
                    prompt = 'Select video data (CleverSys)';
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
        
        function onDoricFilenameEdit(obj)
            target = obj.h.doricFilenameEdit;
            filename = target.Value;
            try
                channels = csvHeader(filename);
                channels{1} = 'time';
                nChannels = numel(channels);
                for i = 1:nChannels
                    channels{i} = sprintf('#%i %s', i, channels{i});
                end
                success = true;
            catch
                success = false;
            end
            setSuccessColor(target, success);
            if success
                obj.doricChannelsEntries = channels;
                updateList(obj.h.doricChannelsList, obj.doricChannelsEntries, obj.settings.doricChannelsChoices);
                obj.saveSettings('doricFilename', filename);
            else
                updateList(obj.h.doricChannelsList, {}, obj.settings.doricChannelsChoices);
            end
        end
        
        function onCleverSysSheetDrop(obj, varargin)
            if numel(varargin) == 1
                sheet = varargin{1};
                epochs = loadCleverSys(obj.h.cleverSysFilenameEdit.Value, sheet);
                entries = epochs(1:2:end);
                obj.cleverSysEpochsEntries = entries;
                updateList(obj.h.cleverSysEpochsList, obj.cleverSysEpochsEntries, obj.settings.cleverSysEpochsChoices);
                obj.saveSettings('cleverSysSheet', sheet);
            else
                updateList(obj.h.doricChannelsList, {}, obj.settings.doricChannelsChoices);
            end
        end
        
        function onCleverSysFilenameEdit(obj)
            target = obj.h.cleverSysFilenameEdit;
            filename = target.Value;
            try
                [~, sheets, tableFormat] = xlsfinfo(filename);
                success = true;
            catch
                success = false;
            end
            setSuccessColor(target, success);
            if success
                obj.cleverSysIsCSV = ismember(tableFormat, {'xlHtml', 'xlCSV'});
                obj.h.cleverSysSheetDrop.Items = sheets;
                if ismember(obj.settings.cleverSysSheet, sheets)
                    sheet = obj.settings.cleverSysSheet;
                else
                    sheet = sheets{1};
                end
                obj.onCleverSysSheetDrop(sheet);
                obj.saveSettings('cleverSysFilename', filename, 'cleverSysSheet', sheet);
            else
                obj.h.cleverSysSheetDrop.Items = {};
                obj.onCleverSysSheetDrop();
            end
        end
        
        function onLabChartFilenameEdit(obj)
            target = obj.h.labChartFilenameEdit;
            filename = target.Value;
            try
                [information, ~, comments] = Aditch.getInformation(filename);
                success = true;
            catch
                success = false;
            end
            setSuccessColor(target, success);
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
                updateList(obj.h.labChartChannelsList, obj.labChartChannelsEntries, obj.settings.labChartChannelsChoices);
                obj.labChartBlocksEntries = blocks;
                updateList(obj.h.labChartBlocksList, obj.labChartBlocksEntries, obj.settings.labChartBlocksChoices);
                obj.saveSettings('labChartFilename', filename);
            else
                obj.h.labChartCommentsHtml.HTMLSource = '<div></div>';
                updateList(obj.h.labChartChannelsList, {}, obj.settings.labChartChannelsChoices);
                updateList(obj.h.labChartBlocksList, {}, obj.settings.labChartBlocksChoices);
            end
        end
        
        function onBleachingEpochsEdit(obj)
            target = obj.h.fpBleachingEpochsEdit;
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
                obj.saveSettings('fpBleachingEpochs', text);
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
                obj.conditionEpochs = epochs;
                obj.saveSettings('conditionEpochs', text);
            end
        end
        
        function setChoice(obj, list, choice)
            namesShowing = list.Items;
            index = ismember(namesShowing, list.Value);
            
            switch list
                case obj.h.doricChannelsList
                    [obj.settings.doricChannelsChoices{index}] = deal(choice);
                    updateList(list, obj.doricChannelsEntries, obj.settings.doricChannelsChoices);
                    obj.saveSettings('doricChannelsChoices', obj.settings.doricChannelsChoices);
                case obj.h.cleverSysEpochsList
                    [obj.settings.cleverSysEpochsChoices{index}] = deal(choice);
                    updateList(list, obj.cleverSysEpochsEntries, obj.settings.cleverSysEpochsChoices);
                    obj.saveSettings('cleverSysEpochsChoices', obj.settings.cleverSysEpochsChoices);
                case obj.h.labChartChannelsList
                    [obj.settings.labChartChannelsChoices{index}] = deal(choice);
                    updateList(list, obj.labChartChannelsEntries, obj.settings.labChartChannelsChoices);
                    obj.saveSettings('labChartChannelsChoices', obj.settings.labChartChannelsChoices);
                case obj.h.labChartBlocksList
                    [obj.settings.labChartBlocksChoices{index}] = deal(choice);
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
end

function settings = defaults()
    % Internal.
    settings.lastFolder = fullfile(getenv('USERPROFILE'), 'Documents');
    
    % General.
    settings.epochsType = 'Manual';
            
    % Frequency to resample EMG and FP.
    settings.resamplingFrequency = 100;
    % Duration of behavior after first camera pulse.
    settings.behaviorDuration = 1800;
    % Maximum cross-correlation lag between EMG and FP.
    settings.xcorrLag = 5;
    
    % EMG.
    % How to filter EMG.
    settings.emgBandpassType = 'Zero-phase';
    % How to calculate envelope.
    settings.emgEnvelopeType = 'Low-pass';
    % Frequency of the envelope filter when choice is low-pass.
    settings.emgEnvelopeLowpassFrequency = 5;
    % Duration of the envelope window when choice is rms.
    settings.emgEnvelopeSize = 0.9;
    % Burden et al 2003. Muscles: vastus lateralis, vastus medialis, semitendinosus, biceps femoris.
    %   Under 20Hz: noise.
    %   Under 450Hz: still clearly identifiable rapid on-off bursts of the EMG.
    settings.emgBandpassFrequencyLow = 20;
    settings.emgBandpassFrequencyHigh = 120;
    
    % FP.
    settings.fpFilename = '';
    settings.fpBleachingEpochs = '[-Inf, Inf]';
    settings.fpPeaksBandpassFrequencyLow = 0.02;
    settings.fpPeaksBandpassFrequencyHigh = 0.2;
    settings.fpDffLowpassFrequency = settings.resamplingFrequency;
    settings.fpBleachingLowpassFrequency = 0.1;
    settings.fpThresholdingFunction = 'mad';
    settings.fpThresholdingFactor = 2;
    settings.fpNormalizationType = 'z-score';
    settings.fpWindowSize = 120;
    
    % Epochs.
    settings.conditionEpochs = '';
    
    % CleverSys epochs.
    settings.cleverSysFilename = '';
    settings.cleverSysSheet = '';
    settings.cleverSysEpochsChoices = {};
    
    % Doric.
    settings.doricFilename = '';
    settings.doricChannelsChoices = {};
    
    % LabChart.
    settings.labChartFilename = '';
    settings.labChartChannelsChoices = {};
    settings.labChartBlocksChoices = {};
    
    % Plots.
    settings.plotFpTrace = true;
    settings.plotFpPower = true;
    settings.plotFpStats = true;
    settings.plotFpTrigger = true;
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
        
function updateList(list, entries, choices)
    labels = cell(size(entries));
    charCounts = cellfun(@numel, entries);
    tabCounts = max(charCounts) - charCounts + 2;
    nEntries = numel(entries);
    for i = 1:nEntries
        entry = entries{i};
        if numel(choices) >= i && ischar(choices{i})
            choice = choices{i};
        else
            choice = 'Unset';
        end
        if isequal(lower(choice), 'unset')
            label = entry;
        else
            label = sprintf('%s%*s==> %s', entry, tabCounts(i), ' ', choice);
        end
        labels{i} = label;
    end
    index = ismember(list.Items, list.Value);
    index(nEntries + 1:end) = [];
    list.Items = labels;
    list.Value = labels(index);
end

function setSuccessColor(target, success)
    if success
        target.BackgroundColor = [1.0, 1.0, 1.0];
    else
        target.BackgroundColor = [1.0, 0.5, 0.5];
    end
end

function header = csvHeader(filename)
    fid = fopen(filename, 'r');
    line = fgetl(fid);
    header = textscan(line, '%s', 'Delimiter', ',');
    header = header{1};
    fclose(fid);
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

%#ok<*PROPLC>