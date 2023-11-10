% data = Doric.load(filename)
% Load Doric data in with h5 or doric extension.

% 2023-02-03. Leonardo Molina.
% 2023-11-09. Last modified.
function data = load(filename, datasets)
    if nargin < 2
        datasets1 = {
            '/DataAcquisition/FPConsole/Signals/Series0001/AIN01xAOUT01-LockIn/Time'
            '/DataAcquisition/FPConsole/Signals/Series0001/AIN01xAOUT01-LockIn/Values'
            '/DataAcquisition/FPConsole/Signals/Series0001/AIN02xAOUT02-LockIn/Values'
        };
        datasets2 = {
            '/Traces/Console/Time(s)/Console_time(s)'
            '/Traces/Console/AIn-1 - Dem (AOut-1)/AIn-1 - Dem (AOut-1)'
            '/Traces/Console/AIn-2 - Dem (AOut-2)/AIn-2 - Dem (AOut-2)'
        };
        try
            h5info(filename, datasets1{1});
            datasets = datasets1;
        catch
            datasets = datasets2;
        end
    end
    nDatasets = numel(datasets);
    for i = 1:nDatasets
        column = h5read(filename, datasets{i});
        if i == 1
            nRows = numel(column);
            data = zeros(nRows, nDatasets);
        end
        data(:, i) = column;
    end
    k = ~any(isnan(data), 2);
    data = data(k, :);
end