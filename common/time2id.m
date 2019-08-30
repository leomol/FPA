% id = time2id(timeVector, epochs)
% Returns the index in timeVector where timeVector is enclosed by the given
% epochs (from1, to1, from2, to2).

% 2019-02-01. Leonardo Molina.
% 2019-08-30. Last modified.
function id = time2id(timeVector, epochs)
    epochs = epochs(:);
    nEpochs = numel(epochs) / 2;
    epochs = zeros(2, nEpochs);
    epochs(1:2:end) = arrayfun(@(l) find(timeVector >= l, 1, 'first'), epochs(1:2:end));
    epochs(2:2:end) = arrayfun(@(l) find(timeVector <= l, 1, 'last'), epochs(2:2:end));
    id = arrayfun(@(e) epochs(1, e):epochs(2, e), 1:nEpochs, 'UniformOutput', false);
    id = [id{:}];
end