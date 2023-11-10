% [data, names] = TDT.load(folder, channels)
% Parses a project folder recorded with TDT DAQ / Synapse and returns a data
% matrix where columns correspond to channels listed in names.

% 2019-02-01. Leonardo Molina.
% 2023-09-11. Last modified.
function data = load(folder, names)
    raw = TDTbin2mat(folder, 'TYPE', {'epocs', 'scalars', 'streams'}, 'CHANNEL', 1);
    frequency = raw.streams.(names{1}).fs;
    nSamples = numel(raw.streams.(names{1}).data);
    nChannels = numel(names);
    data = NaN(nSamples, nChannels + 1);
    for c = 1:nChannels
        data(:, c + 1) = double(raw.streams.(names{c}).data);
    end
    data(:, 1) = transpose(1:nSamples) / frequency;
end