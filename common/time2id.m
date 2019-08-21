function id = time2id(timeVector, timeLimits)
    timeLimits = timeLimits(:);
    nEpochs = numel(timeLimits) / 2;
    epochs = zeros(2, nEpochs);
    epochs(1:2:end) = arrayfun(@(l) find(timeVector >= l, 1, 'first'), timeLimits(1:2:end));
    epochs(2:2:end) = arrayfun(@(l) find(timeVector <= l, 1, 'last'), timeLimits(2:2:end));
    id = arrayfun(@(e) epochs(1, e):epochs(2, e), 1:nEpochs, 'UniformOutput', false);
    id = [id{:}];
end