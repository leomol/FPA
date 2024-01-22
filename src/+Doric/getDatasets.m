% datasets = Doric.getDatasets(filename)
% List datasets within Doric's h5 datafile.

% 2023-02-03. Leonardo Molina.
% 2024-01-22. Last modified.
function datasets = getDatasets(filename)
    root = h5info(filename, '/');
    datasets = get(filename, root);
    datasets = datasets.cellstr();
end

function locations = get(filename, location)
    nGroups = numel(location.Groups);
    nDatasets = numel(location.Datasets);
    locations = cell(0, 1);
    if nGroups > 0
        for i = 1:nGroups
            locations = cat(1, locations, get(filename, location.Groups(i)));
        end
    elseif nDatasets > 0
        locations = location.Name + "/" + {location.Datasets.Name}';
    end
end