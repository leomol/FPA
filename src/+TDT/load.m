% [data, names] = TDT.load(folder, channels)
% Parses a project folder recorded with TDT DAQ / Synapse and returns a data
% matrix where columns correspond to channels listed in names.
% Includes all streams listed as long as they match in size with the first stream.

% 2019-02-01. Leonardo Molina.
% 2025-11-24. Last modified.
function [data, names] = load(folder, names)
    raw = TDTbin2mat(folder, 'TYPE', {'epocs', 'scalars', 'streams'}, 'CHANNEL', 1);
    if nargin < 2
        names = fieldnames(raw.streams)';
    end
    frequency = raw.streams.(names{1}).fs;
    nSamples = numel(raw.streams.(names{1}).data);
    nChannels = numel(names);
    data = NaN(nSamples, nChannels + 1);
    mask = true(1, nChannels);
    for c = 1:nChannels
        d = double(raw.streams.(names{c}).data);
        if numel(d) == nSamples
            data(:, c + 1) = d;
        else
            mask(c) = false;
        end
    end
    data(:, 1) = transpose(1:nSamples) / frequency;
    data = data(:, [true, mask]);
    names = ['time', names(mask)];
end