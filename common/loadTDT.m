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