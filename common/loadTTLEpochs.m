% [lowEpochs, highEpochs] = loadTTLEpochs(filename, column)
% For a recording where TTL=0 encodes condition 1 and TTL=1 encodes condition 2
% lowEpochs and highEpochs are timestamp pairs separated for each condition.
% 
% 2020-02-03. Leonardo Molina.
% 2020-02-10. Last modified.
function [lowEpochs, highEpochs] = loadTTLEpochs(filename, column)
    data = loadData(filename);
    time = data(:, 1);
    ttl = data(:, column) > 0;
    delta = diff(ttl);
    % Low starts at t=-Inf, continues until rise, then restarts at the next fall.
    lowEpochs = time([delta == +1; false] | [false; delta == -1]);
    highEpochs = time([false; delta == +1] | [delta == -1; false]);
    % Complete epochs according to last change.
    if mod(numel(lowEpochs) + 1, 2) == 0
        highEpochs(end + 1) = Inf;
    else
        lowEpochs(end + 1) = Inf;
    end
    lowEpochs = [-Inf; lowEpochs];
    nLow = numel(lowEpochs) / 2;
    nHigh = numel(highEpochs) / 2;
    lowEpochs = reshape(lowEpochs, 2, nLow);
    highEpochs = reshape(highEpochs, 2, nHigh);
end