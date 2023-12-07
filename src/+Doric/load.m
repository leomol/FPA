% data = Doric.load(filename)
% Load Doric data in with h5 or doric extension using defaults.
% 
% patterns = {'.*Console_time.*', '.*01-LockIn/Values', '.*02-LockIn/Values'};
% data = Doric.load(filename, patterns)
% Repeat but read data from the given locations.

% 2023-02-03. Leonardo Molina.
% 2023-12-07. Last modified.
function data = load(filename, patterns)
    if nargin < 2
        patterns = {
            % /DataAcquisition/FPConsole/Signals/Series0001/AIN01xAOUT01-LockIn/Time
            % /DataAcquisition/EPConsole/Signals/Series0001/AIN01xAOUT02-LockIn/Time
            % /DataAcquisition/EPConsole/Signals/Series0001/AnalogIn/Time
            % /Traces/Console/Time(s)/Console_time(s)
            '.*(LockIn/Time|/Time|time\(s\))$'

            % /DataAcquisition/EPConsole/Signals/Series0001/AIN01xAOUT01-LockIn/Values
            % /DataAcquisition/EPConsole/Signals/Series0001/AnalogIn/AIN01
            % /Traces/Console/AIn-1 - Dem (AOut-1)/AIn-1 - Dem (AOut-1)
            '.*/(AIN01xAOUT01-LockIn/Values|AnalogIn/AIN01|AIn-1 - Dem.*)$'

            % /DataAcquisition/FPConsole/Signals/Series0001/AIN01xAOUT02-LockIn/Values
            % /DataAcquisition/EPConsole/Signals/Series0001/AnalogIn/AIN02
            % /Traces/Console/AIn-2 - Dem (AOut-2)/AIn-2 - Dem (AOut-2)
            '.*/(AIN01xAOUT02-LockIn/Values|AnalogIn/AIN02|AIn-2 - Dem.*)$'
            };
    end
    nTargets = numel(patterns);
    info = h5info(filename, '/DataAcquisition');
    locations = getLocations(filename, info);
    match = @(pattern) locations(find(~cellfun(@isempty, regexp(locations, pattern, 'match')), 1));
    data = h5read(filename, match(patterns{1}));
    data(:, 2:nTargets) = NaN;
    for i = 2:nTargets
        location = match(patterns{i});
        if isempty(location)
            error('[Doric.read] Could not match pattern %i against available locations:\n%s', i, sprintf('  %s\n', locations{:}));
        else
            data(:, i) = h5read(filename, match(patterns{i}));
        end
    end
end

function locations = getLocations(filename, dataAcquisition)
    nGroups = numel(dataAcquisition.Groups);
    nDatasets = numel(dataAcquisition.Datasets);
    if nGroups > 0
        locations = zeros(0, 1);
        for i = 1:nGroups
            locations = cat(1, locations, getLocations(filename, dataAcquisition.Groups(i)));
        end
    elseif nDatasets > 0
        locations = dataAcquisition.Name + "/" + string({dataAcquisition.Datasets.Name}');
    end
end