% id = time2id(timeVector, epochs)
% Returns the index in timeVector where timeVector is enclosed by the given
% epochs (from1, to1, from2, to2).

% 2019-02-01. Leonardo Molina.
% 2019-10-04. Last modified.
function id = time2id(timeVector, epochs)
    epochs = epochs(:);
    nEpochs = numel(epochs) / 2;
    timeLimits = zeros(2, nEpochs);
    timeLimits(1:2:end) = arrayfun(@(l) find(timeVector >= l, 1, 'first'), epochs(1:2:end));
    timeLimits(2:2:end) = arrayfun(@(l) find(timeVector <= l, 1, 'last'), epochs(2:2:end));
    id = arrayfun(@(e) colon(timeLimits(1, e),timeLimits(2, e))', 1:nEpochs, 'UniformOutput', false);
    id = cat(1, id{:});
end