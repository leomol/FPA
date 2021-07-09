% epochs = loadBinaryStates(filename)
% filename is a csv file where columns are time, frame number, and one extra column for each behavior.
% Each behavior column contains either 0's or 1's, encoding presence or absence of such behavior for each frame.
% Special cases:
%   Multiple behaviors may occur at once.
%   A behavior may last a single frame.
%   Blanks are re-encoded as 0's.
% 
% For example:
%   time, frame, approach, rearing
%   0.00,     1,        0,       0
%   0.20,     2,        1,       0
%   0.40,     3,        1,       0
%   0.60,     4,        1,       1
%   0.80,     5,        1,       1
%   1.00,     6,        0,       1
%   1.20,     7,        0,       1
%   1.40,     8,        1,       0
%   1.60,     9,        0,       0
% 
% {'approach', [0.20, 0.80, 1.40, 1.40], 'rearing', [0.60, 1.20]}

% 2021-07-08. Leonardo Molina.
% 2021-07-08. Last modified.
function epochs = loadBinaryStates(filename)
    data = readtable(filename);
    % Retrieve behavioral labels. First two are time and frame number.
    labels = data.Properties.VariableNames(3:end);
    % Turn table into matrix.
    time = data.time;
    data = data{:, 3:end};
    % Re-encode blanks to zeros.
    data(isnan(data)) = 0;
    nBehaviors = size(data, 2);
    % Initialize epochs cell.
    epochs = cell(1, 2 * nBehaviors);
    epochs(1:2:end) = labels;
    for c = 1:nBehaviors
        % For each behavior column, find changes.
        b = data(:, c);
        d = diff(b) ~= 0;
        f = find(d);
        % Behavior starts with the change.
        f(1:2:end) = f(1:2:end) + 1;
        % Return time of start/end of the change.
        epochs{2 * c} = time(f);
    end
end
