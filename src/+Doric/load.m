% data = Doric.load(filename)
% Load Doric's hdf5 datafiles (.doric) using typical dataset locations
% corresponding to time, and demodulated signals 1 and 2.
% 
% datasets = {'01-LockIn/Time', '01-LockIn/Values', '02-LockIn/Values'};
% data = Doric.load(filename, datasets)
% Same as before but provide the name of the datasets (interpreted as regular expressions).
% 
% datasets = {
%   '/Traces/Console/Time(s)/Console_time(s)'
%   '/Traces/Console/AIn-1 - Dem (AOut-1)/AIn-1 - Dem (AOut-1)'
%   '/Traces/Console/AIn-2 - Dem (AOut-2)/AIn-2 - Dem (AOut-2)'
% };
% data = Doric.load(filename, datasets, false)
% Same as before but the name of the datasets are interpreted literally (as opposed to regular expressions).

% 2023-02-03. Leonardo Molina.
% 2024-01-19. Last modified.
function data = load(filename, datasets, regex)
    if nargin < 2
        datasets = {
            % AIN01xAOUT01-LockIn/Time
            % AIN02xAOUT02-LockIn/Time
            % AnalogIn/Time
            % Console_time(s)
            '.*(LockIn/Time|AnalogIn/Time|Console_time\(s\))$'
            
            % AIN01xAOUT01-LockIn/Values
            % AnalogIn/AIN01
            % AIn-1 - Dem (AOut-1)
            '.*/(AIN01xAOUT01-LockIn/Values|AnalogIn/AIN01|AIn-1 - Dem.*)$'
            
            % AIN02xAOUT02-LockIn/Values
            % AnalogIn/AIN02
            % AIn-2 - Dem (AOut-2)
            '.*/(AIN02xAOUT02-LockIn/Values|AnalogIn/AIN02|AIn-2 - Dem.*)$'
            };
    end
    if nargin < 3
        regex = true;
    end
    
    nTargets = numel(datasets);
    available = Doric.getDatasets(filename);

    unavailable = {};
    targets = cell(size(datasets));
    
    if regex
        for i = 1:nTargets
            matches = regexp(available, datasets{i}, 'match');
            found = ~cellfun(@isempty, matches);
            if any(found)
                k = find(found, 1);
                targets{i} = available{k};
            else
                unavailable = cat(1, unavailable, datasets{i});
            end
        end
    else
        for i = 1:nTargets
            if ismember(datasets{i}, available)
                targets{i} = datasets{i};
            else
                unavailable = cat(1, unavailable, datasets{i});
            end
        end
    end

    if numel(unavailable) > 0
        if regex
            msg = 'Could not match the following patterns against any available locations in the dataset';
        else
            msg = 'The dataset does not contain the locations';
        end
        msg = sprintf('%s:\n  "%s"', msg, strjoin(unavailable, '"\n  "'));
        msg = sprintf('%s\n\nThe dataset only has the following locations available:\n  "%s"', msg, strjoin(available, '"\n  "'));
        error(msg, '');
    end
    
    info = h5info(filename, targets{1});
    data = NaN(info.Dataspace.Size, nTargets);
    for i = 1:nTargets
        data(:, i) = h5read(filename, targets{i});
    end
end