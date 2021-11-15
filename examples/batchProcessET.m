% FPA batch processing script for ET's thesis data.
% This script will load Doric csv files and their corresponding behavioral data (Murari Lab spec), all of which are provided in the project variable below.
% Data is processed according to user-defined settings, then triggered data is generated in the form of output csv files and plots.
% 
% Note:
%   -Do not use dots in folder or file names (IMPORTANT).
%   -The configuration was defined according to a conversation in Slack (Chat Elizabeth+Patrick+Leo 2021-11-14).
%   -Homework:
%     -Check all settings (including data columns), folders and epoch defitions below!
%     -Add new data (if any) following the pattern below. Ask if unsure.
% 
% 2021-09-30. Leonardo Molina.
% 2021-11-15. Last modified.

% Top project folder.
top = 'C:\Users\Molina\Documents\public\data\HALO\FibrePhotometry\Elizabeth\Thesis';
% Formats to exports images to. A single or multiple image formats.
formats = {'png', 'fig', 'epsc'};
% Suffix added to the exported folder.
outputSuffix = ' - exported';

% General configuration to apply to all files.
configuration = struct();
configuration.resamplingFrequency = 100;
configuration.baselineLowpassFrequency = 0.1;
configuration.fitReference = false;
configuration.lowpassFrequency = 4;
configuration.peaksLowpassFrequency = 0.001;
configuration.threshold = {2.91, @mad, @median};

% Option 1: z-score.
configuration.f0 = @mean;
configuration.f1 = @std;

% Option 2: 5min moving window using median and mad.
% configuration.f0 = {@median, 300};
% configuration.f1 = {@mad, 300};

% Columns corresponding to 465nm and 405nm.
signalColumn = 5;
referenceColumn = 3;

project = { ...
		"familiar male", ...
		[...
            dir(fullfile(top, "Experiment 1.5*\familiar male*\*.csv"))
            dir(fullfile(top, "Experiment 1.6*\familiar male*\*.csv"))
            dir(fullfile(top, "Experiment 1.7*\familiar male*\*.csv"))
        ], ...
        [...
            dir(fullfile(top, "Experiment 1.5*\familiar male*\tracking*\**\*.csv"))
            dir(fullfile(top, "Experiment 1.6*\familiar male*\tracking*\**\*.csv"))
            dir(fullfile(top, "Experiment 1.7*\familiar male*\tracking*\**\*.csv"))
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
		"novel female", ...
		[...
            dir(fullfile(top, "Experiment 1.5*\novel female*\*.csv"))
            dir(fullfile(top, "Experiment 1.6*\novel female*\*.csv"))
            dir(fullfile(top, "Experiment 1.7*\novel female*\*.csv"))
        ], ...
        [...
            dir(fullfile(top, "Experiment 1.5*\novel female*\tracking*\**\*.csv"))
            dir(fullfile(top, "Experiment 1.6*\novel female*\tracking*\**\*.csv"))
            dir(fullfile(top, "Experiment 1.7*\novel female*\tracking*\**\*.csv"))
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
		"novel male", ...
		[...
            dir(fullfile(top, "Experiment 1.5*\novel male*\*.csv"))
            dir(fullfile(top, "Experiment 1.6*\novel male*\*.csv"))
            dir(fullfile(top, "Experiment 1.7*\novel male*\*.csv"))
        ], ...
        [...
            dir(fullfile(top, "Experiment 1.5*\novel male*\tracking*\**\*.csv"))
            dir(fullfile(top, "Experiment 1.6*\novel male*\tracking*\**\*.csv"))
            dir(fullfile(top, "Experiment 1.7*\novel male*\tracking*\**\*.csv"))
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
		"novel object", ...
		[...
            dir(fullfile(top, "Experiment 1.5*\novel object*\*.csv"))
            dir(fullfile(top, "Experiment 1.6*\novel object*\*.csv"))
            dir(fullfile(top, "Experiment 1.7*\novel object*\*.csv"))
        ], ...
        [...
            dir(fullfile(top, "Experiment 1.5*\novel object*\tracking*\**\*.csv"))
            dir(fullfile(top, "Experiment 1.6*\novel object*\tracking*\**\*.csv"))
            dir(fullfile(top, "Experiment 1.7*\novel object*\tracking*\**\*.csv"))
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
            config.conditionEpochs = loadBinaryStates(inputEventFile);

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