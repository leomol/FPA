% ttl = loadInscopixTTL(filename, ttlTarget)
% Returns timestamps where tllTarget changes from low to high and high to low in Inscopix DAQ.

% 2019-02-01. Leonardo Molina.
% 2020-02-19. Last modified.
function [rise, fall] = loadInscopixTTL(filename, ttlTarget)
    % Get time of transitions from low to high.
    fid = fopen(filename, 'r');
    parts = textscan(fid, '%f%s%d', 'Delimiter', ',', 'HeaderLines', 1);
    fclose(fid);
    if nargin < 2
        ttlTarget = parts{2}{1};
    end
    k = ismember(parts{2}, ttlTarget);
    time = parts{1}(k);
    state = parts{3}(k);
    rise = time([false; diff(state) == +1]);
    fall = time([false; diff(state) == -1]);
    rise = rise(:);
    fall = fall(:);
end