% [faces, vertices] = patchEpochs(epochs, minimum, maximum)
% Create vectors faces and vertices from vector epochs (from1, to1, from2, to2, ...)
% to feed to MATLAB's function patch. The y-limits of the "drawing" are
% given by minimum and maximum.

% 2019-02-01. Leonardo Molina.
% 2019-08-30. Last modified.
function [faces, vertices] = patchEpochs(epochs, minimum, maximum)
    epochs = epochs(:);
    nEpochs = numel(epochs);
    vertices = zeros(2, 2 * nEpochs);
    for e = 1:2:nEpochs
        range = 4 * (e - 1) + (1:8);
        e1 = e;
        e2 = e + 1;
        vertices(range) = [epochs(e1), minimum, epochs(e1), maximum, epochs(e2), maximum, epochs(e2), minimum];
    end
    vertices = vertices';
    faces = reshape(1:2 * nEpochs, 4, nEpochs / 2)';
end