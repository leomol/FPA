function [faces, vertices] = patchEpochs(epochs, mn, mx)
    epochs = epochs(:);
    nEpochs = numel(epochs);
    vertices = zeros(2, 2 * nEpochs);
    for e = 1:2:nEpochs
        range = 4 * (e - 1) + (1:8);
        e1 = e;
        e2 = e + 1;
        vertices(range) = [epochs(e1), mn, epochs(e1), mx, epochs(e2), mx, epochs(e2), mn];
    end
    vertices = vertices';
    faces = reshape(1:2 * nEpochs, 4, nEpochs / 2)';
end