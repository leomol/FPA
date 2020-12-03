% id = time2id(timeVector, epochs)
% Returns the index in timeVector where timeVector is enclosed by the given
% epochs (from1, to1, from2, to2).

% 2019-02-01. Leonardo Molina.
% 2020-11-24. Last modified.
function id = time2id(timeVector, epochs)
    epochs = epochs(:);
    % timeLimits = zeros(2, nEpochs);
    a = arrayfun(@(k) find(timeVector >= k, 1, 'first'), epochs(1:2:end), 'UniformOutput', false);
    b = arrayfun(@(k) find(timeVector <= k, 1, 'last'), epochs(2:2:end), 'UniformOutput', false);
    k = cellfun(@isempty, a) | cellfun(@isempty, b);
    a(k) = [];
    b(k) = [];
    nEpochs = numel(a);
    timeLimits = zeros(2, nEpochs);
    timeLimits(1:2:end) = [a{:}];
    timeLimits(2:2:end) = [b{:}];
    id = arrayfun(@(e) colon(timeLimits(1, e), timeLimits(2, e))', 1:nEpochs, 'UniformOutput', false);
    id = cat(1, id{:});
end