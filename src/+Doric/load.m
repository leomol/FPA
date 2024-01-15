% data = Doric.load(filename)
% Load Doric data in with h5 or doric extension using default locations.
% 
% datasets = {'.*01-LockIn/Time', '.*01-LockIn/Values', '.*02-LockIn/Values'};
% data = Doric.load(filename, datasets)
% Repeat but read data from the given datasets.

% 2023-02-03. Leonardo Molina.
% 2024-01-18. Last modified.
function data = load(filename, datasets)
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
    nTargets = numel(datasets);
    locations = Doric.getDatasets(filename);
    match = @(pattern) locations(find(~cellfun(@isempty, regexp(locations, pattern, 'match')), 1));
    data = h5read(filename, match(datasets{1}));
    data(:, 2:nTargets) = NaN;
    for i = 2:nTargets
        location = match(datasets{i});
        if isempty(location)
            error('[Doric.load] Could not match pattern %i against available locations:\n%s', i, sprintf('  %s\n', locations{:}));
        else
            data(:, i) = h5read(filename, location);
        end
    end
end