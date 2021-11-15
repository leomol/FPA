% batchProcessCB(projectFolder, <dataFileFilter>, <eventFileFilter>, <configuration>, <override>, <formats>, <suffix>)
% Scans projectFolder recursively, creates an output folder matching FP files, exports FPA data and figures to this folder.
% 
%   projectFolder: root folder with arbitrary organization but where FP files
%      and event files must reside within the same folder.
%   dataFileFilter: filter FP files by name. For example:
%     '*.csv'
%     '*Excursion.csv'
%     ...
%   eventFileFilter: filter event files by name. For example:
%     '*.tsv'
%   configuration: FPA configuration structure.
%   override: Whether to override previously exported data.
%   formats: cell array with image types to export. For example:
%     {'png'}
%     {'png', 'fig', 'epsc'}
%   suffix: append this string to the output folder for revision control.
% 
% Example 1 - use all defaults:
%   projectFolder = 'C:\Users\Molina\OneDrive - University of Calgary\Sandeep FP Experiment';
%   batchProcessCB(projectFolder);
% 
% Example 2:
%   projectFolder = 'C:\Users\Molina\OneDrive - University of Calgary\Sandeep FP Experiment';  
%   dataFileFilter = '*Excursion.csv';
%   eventFileFilter = '*.tsv';
%   override = true;
%   formats = {'png', 'fig', 'epsc'};
%   suffix = ' - exported';
% 
%   configuration = struct();
%   configuration.lowpassFrequency = 2;
%   configuration.f0 = @median;
%   configuration.f1 = @mad;
%   configuration.threshold = {2.91, @mad, @median};
%   configuration.baselineEpochs = [-Inf, Inf];
%   
%   batchProcessCB(projectFolder, dataFileFilter, eventFileFilter, configuration, override, formats, suffix);
% 
% Note:
%   -This script currently assumes that data files were acquired with Doric DAQ
%    and that event files were annotated with BORIS.
%   -Do not use dots in folder or file names.

% 2021-09-16. Leonardo Molina.
% 2021-09-16. Last modified.
function batchProcessCB(projectFolder, dataFileFilter, eventFileFilter, configuration, override, formats, outputSuffix)
    % Add dependencies.
    folder = fileparts(mfilename('fullpath'));
    addpath(fullfile(folder, '..'));
    addpath(genpath(fullfile(folder, '../common')));
    
    if nargin < 2
        dataFileFilter = '*.csv';
    end
    
    if nargin < 3
        eventFileFilter = '*.tsv';
    end
    
    if nargin < 4
        configuration = struct();
    end
    
    if nargin < 5
        override = false;
    end
    
    if nargin < 6
        formats = {'png'};
    end
    
    if nargin < 7
        outputSuffix = '';
    end
    
    % Columns corresponding to 465nm and 405nm.
    signalColumn = 5;
    referenceColumn = 3;

    % Retrieve all files.
    dataFiles = dir(fullfile(projectFolder, '**', dataFileFilter));
    eventFiles = dir(fullfile(projectFolder, '**', eventFileFilter));
    
    % Only process matching files on each folder.
    dataUID = {dataFiles.folder};
    eventsUID = {eventFiles.folder};
    [~, k1, k2] = intersect(dataUID, eventsUID);
    dataFiles = dataFiles(k1);
    eventFiles = eventFiles(k2);
    dataFiles = fullfile({dataFiles.folder}, {dataFiles.name});
    eventFiles = fullfile({eventFiles.folder}, {eventFiles.name});
    nFiles = numel(dataFiles);
    
    % Iterate.
    userFigures = findobj('Type', 'figure');
    for i = 1:nFiles
        inputDataFile = dataFiles{i};
        inputEventFile = eventFiles{i};
        
        % Create output folder.
        [folder, basename] = fileparts(inputDataFile);
        outputFolder = fullfile(folder, [basename, outputSuffix]);
        do = true;
        fprintf('\n[%i:%i] %s\n', i, nFiles, inputDataFile);
        if exist(outputFolder, 'dir') == 7
            if override
                fprintf('Overriding data.\n');
            else
                do = false;
                fprintf('Output folder already exists.\n');
            end
        else
            mkdir(outputFolder);
        end

        if do
            % Fiber-photometry data recorded with Doric DAQ.
            data = loadData(inputDataFile);
            time = data(:, 1);
            signal = data(:, signalColumn);
            reference = data(:, referenceColumn);

            % Event files scored with Boris.
            configuration.conditionEpochs = loadBoris(inputEventFile);

            % Call FPA with given configuration.
            fpa = FPA(time, signal, reference, configuration);
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
    end
end