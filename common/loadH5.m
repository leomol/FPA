% data = loadDoric(filename)
% Load Doric data in with h5 or doric extension.

% 2023-02-03. Leonardo Molina.
% 2023-02-03. Last modified.
function data = loadH5(filename, datasets)
    if nargin < 2
        datasets = {
            '/DataAcquisition/FPConsole/Signals/Series0001/AIN01xAOUT01-LockIn/Time'
            '/DataAcquisition/FPConsole/Signals/Series0001/AIN01xAOUT01-LockIn/Values'
            '/DataAcquisition/FPConsole/Signals/Series0001/AIN02xAOUT02-LockIn/Values'
            };
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
end