% [rise, fall] = loadTTL(filename, <column>)
% Returns timestamps where pin changes from low to high and high to low.
%
% 2020-02-03. Leonardo Molina.
% 2020-01-31. Last modified.
function [rise, fall] = loadTTL(filename, column)
    data = loadData(filename);
    time = data(:, 1);
    ttl = data(:, column);
    delta = diff(ttl);
    rise = time([false; delta == +1]);
    fall = time([false; delta == -1]);
end