% [data, names] = loadTDT(folder)
% Parses a project folder stored by a TDT DAQ and returns a data matrix where
% columns correspond to channels listed in names.

% 2019-02-01. Leonardo Molina.
% 2019-08-30. Last modified.
function [data, names] = loadTDT(folder)
    data = TDTbin2mat(folder, 'TYPE', {'epocs', 'scalars', 'streams'}, 'CHANNEL', 1);
    names = fieldnames(data.streams);
    frequency = data.streams.(names{1}).fs;
    nSamples = numel(data.streams.(names{1}));
    nChannels = numel(names);
    data = NaN(nSamples, nChannels + 1);
    for c = 1:nChannels
        data(:, c + 1) = double(data.streams.(names{c}).data);
    end
    data(:, 1) = step * transpose(1:nSamples) / frequency;
    names = ['names'; names];
end