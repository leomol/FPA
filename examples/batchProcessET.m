% FPA batch processing script for ET's thesis data.
% This script will load Doric csv files and their corresponding behavioral data (CleverSys or BinaryStates), all of which are provided in the project variable below.
% Data is processed according to user-defined settings, then triggered data is generated in the form of output csv files and plots.
% 
% Note:
%   -Do not use dots in folder or file names (IMPORTANT).
%   -Add new data (if any) by following the pattern in the project variable.
%   -Configuration and processing steps have been re-defined based on on-going conversations:
%      2021-11-14 Slack chat with Elizabeth and Patrick.
%      2021-11-15 Anydesk messaging with Elizabeth: peaksLowpassFrequency was reverted back to the defaults (0.5Hz).
%      2021-11-18 Slack: Load both CleverSys or BinaryStates.
%   -Always check:
%      All settings including data columns, folders and epoch defitions.
% 
% 2021-09-30. Leonardo Molina.
% 2021-11-24. Last modified.

% Top project folder.
top = 'G:\My Drive\MSc\Research Data\Fiber Photometry Experiments\';
% Formats to exports images to. A single or multiple image formats.
formats = {'png', 'fig'}; % If you wanted eps, add 'epsc'
% Suffix added to the exported folder.
outputSuffix = ' - exported';

% General configuration to apply to all files.
configuration = struct();
configuration.resamplingFrequency = 100;
configuration.baselineLowpassFrequency = 0.1;
configuration.fitReference = false;
configuration.lowpassFrequency = 4;
configuration.peaksLowpassFrequency = 0.5;
configuration.threshold = {2.91, @mad, @median};

% z-score.
configuration.f0 = @mean;
configuration.f1 = @std;

% Columns corresponding to 465nm and 405nm.
signalColumn = 2;
referenceColumn = 4;

search = @(folder) mdir(fullfile(top, folder), {'*.csv', '*.xlsx'});

project = { ...
		'familiar male', ...
		[...
            search('Experiment_1_5*\familiar male*')
            search('Experiment_1_6*\familiar male*')
            search('Experiment_1_7*\familiar male*')
        ], ...
        [...
            search('Experiment_1_5*\familiar male*\tracking*\**')
            search('Experiment_1_6*\familiar male*\tracking*\**')
            search('Experiment_1_7*\familiar male*\tracking*\**')
        ], ...
        [
            struct('baselineEpochs', [ 0, 300, 637, 937], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 655, 955], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 640, 940], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 645, 945], 'artifactEpochs', [])
            struct('baselineEpochs', [   -Inf, Inf     ], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 660, 960], 'artifactEpochs', [654, 670])
            struct('baselineEpochs', [   -Inf, Inf     ], 'artifactEpochs', [])
            ...
            struct('baselineEpochs', [ 0, 300, 635, 935], 'artifactEpochs', [795, 815])
            struct('baselineEpochs', [ 0, 300, 617, 917], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 622, 922], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 617, 917], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 622, 922], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 620, 920], 'artifactEpochs', [])
            ...
            struct('baselineEpochs', [ 0, 300, 650, 950], 'artifactEpochs', [635, 648])
            struct('baselineEpochs', [ 0, 300, 645, 945], 'artifactEpochs', [625, 640])
            struct('baselineEpochs', [ 0, 300, 615, 915], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 617, 917], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 617, 917], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 620, 920], 'artifactEpochs', [])
        ], ...
		...
		'novel female', ...
		[...
            search('Experiment_1_5*\novel female*')
            search('Experiment_1_6*\novel female*')
            search('Experiment_1_7*\novel female*')
        ], ...
        [...
            search('Experiment_1_5*\novel female*\tracking*\**')
            search('Experiment_1_6*\novel female*\tracking*\**')
            search('Experiment_1_7*\novel female*\tracking*\**')
        ], ...
		[
            struct('baselineEpochs', [ 0, 320, 640, 940], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 320, 642, 942], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 320, 635, 935], 'artifactEpochs', [])
            struct('baselineEpochs', [ 7, 307, 640, 940], 'artifactEpochs', [])
            struct('baselineEpochs', [10, 310, 640, 940], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 655, 955], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 632, 932], 'artifactEpochs', [])
            ...
            struct('baselineEpochs', [ 0, 300, 625, 925], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 622, 922], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 615, 915], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 615, 915], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 622, 922], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 622, 922], 'artifactEpochs', [])
            ...
            struct('baselineEpochs', [ 0, 300, 622, 922], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 617, 917], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 617, 917], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 617, 917], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 612, 912], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 625, 925], 'artifactEpochs', [])
		], ...
		...
		'novel male', ...
		[...
            search('Experiment_1_5*\novel male*')
            search('Experiment_1_6*\novel male*')
            search('Experiment_1_7*\novel male*')
        ], ...
        [...
            search('Experiment_1_5*\novel male*\tracking*\**')
            search('Experiment_1_6*\novel male*\tracking*\**')
            search('Experiment_1_7*\novel male*\tracking*\**')
        ], ...
		[
            struct('baselineEpochs', [ 0, 300, 657, 957], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 640, 940], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 650, 950], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 647, 947], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 642, 942], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 675, 975], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 637, 937], 'artifactEpochs', [])
            ...
            struct('baselineEpochs', [ 0, 300, 625, 925], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 622, 922], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 615, 915], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 615, 915], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 622, 922], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 622, 922], 'artifactEpochs', [])
            ...
            struct('baselineEpochs', [ 0, 300, 622, 922], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 617, 917], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 617, 917], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 617, 917], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 612, 912], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 625, 925], 'artifactEpochs', [])
		], ...
		...
		'novel object', ...
		[...
            search('Experiment_1_5*\novel object*')
            search('Experiment_1_6*\novel object*')
            search('Experiment_1_7*\novel object*')
        ], ...
        [...
            search('Experiment_1_5*\novel object*\tracking*\**')
            search('Experiment_1_6*\novel object*\tracking*\**')
            search('Experiment_1_7*\novel object*\tracking*\**')
        ], ...
		[
            struct('baselineEpochs', [ 0, 300, 635,  935], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 635,  935], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 632,  932], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 640,  940], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 635,  935], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 627,  927], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 627,  927], 'artifactEpochs', [])
            ...
            struct('baselineEpochs', [ 0, 300, 615,  915], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 607,  907], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 617,  917], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 607,  907], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 625,  925], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 700, 1000], 'artifactEpochs', [])
            ...
            struct('baselineEpochs', [ 0, 300, 612,  912], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 607,  907], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 615,  915], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 640,  940], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 607,  907], 'artifactEpochs', [])
            struct('baselineEpochs', [ 0, 300, 610,  910], 'artifactEpochs', [])
		]
    };

% Add dependencies.
folder = fileparts(mfilename('fullpath'));
addpath(fullfile(folder, '..'));
addpath(genpath(fullfile(folder, '../common')));

warningType = 'MATLAB:table:ModifiedAndSavedVarnames';
status = warning('query', warningType);
warning('off', warningType);
userFigures = findobj('Type', 'figure');

nConditions = numel(project) / 4;
for c = 1:nConditions
    conditionName = project{4 * c - 3};
    fpFiles = fullfile({project{4 * c - 2}.folder}, {project{4 * c - 2}.name});
    behaviorFiles = fullfile({project{4 * c - 1}.folder}, {project{4 * c - 1}.name});
    configurations = project{4 * c - 0};
    fprintf('\n%02i:%02i "%s"\n', c, nConditions, conditionName);
    
    nFiles = numel(fpFiles);
    if nFiles == numel(behaviorFiles) && nFiles == numel(configurations)
        for f = 1:nFiles
            inputDataFile = fpFiles{f};
            inputEventFile = behaviorFiles{f};
            
            % Create output folder.
            [folder, basename] = fileparts(inputDataFile);
            outputFolder = fullfile(folder, [basename, outputSuffix]);
            fprintf('  [%02i:%02i] %s\n', f, nFiles, inputDataFile);
            if exist(outputFolder, 'dir') ~= 7
                mkdir(outputFolder);
            end
            
            % Merge base configuration with file-specific configuration.
            config = configuration;
            fnames = fieldnames(configurations(f));
            for s = 1:numel(fnames)
                fname = fnames{s};
                config.(fname) = configurations(f).(fname);
            end
            
            % Fiber-photometry data recorded with Doric DAQ.
            data = loadData(inputDataFile);
            time = data(:, 1);
            signal = data(:, signalColumn);
            reference = data(:, referenceColumn);

            % Event files.
            [~, ~, ext] = fileparts(inputEventFile);
            switch lower(ext)
                case 'csv'
                    config.conditionEpochs = loadBinaryStates(inputEventFile);
                case 'xlsx'
                    config.conditionEpochs = loadCleverSys(inputEventFile);
            end

            % Call FPA with given configuration.
            fpa = FPA(time, signal, reference, config);
            % Print warnings.
            if numel(fpa.warnings) > 0
                fprintf(2, 'Warning(s) for "%s":\n', inputDataFile);
                cellfun(@warning, fpa.warnings);
            end

            % Export.
            prefix = fullfile(outputFolder, basename);
            fpa.export(prefix);

            % Plot, save and close.
            fpa.plot();
            allFigures = findobj('Type', 'figure');
            fpaFigures = setdiff(allFigures, userFigures);
            for j = 1:numel(fpaFigures)
                fig = fpaFigures(j);
                for k = 1:numel(formats)
                    cleanName = regexprep(fig.Name, '[^a-zA-Z0-9 ]', '');
                    saveas(fig, fullfile(outputFolder, cleanName), formats{k});
                end
            end
            close(fpaFigures);

            % Save configuration file.
            save(fullfile(outputFolder, 'configuration.mat'), 'configuration');
        end
    else
        error('Mismatch in the number of data files and/or configuration provided.');
    end
end
warning(status.state, warningType);

function list = mdir(folder, filenames)
    list = cellfun(@(filename) dir(fullfile(folder, filename)), filenames, 'UniformOutput', false);
    list = vertcat(list{:});
end