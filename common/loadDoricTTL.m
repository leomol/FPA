% ttl = loadDoricTTL(filename)
% Returns timestamps where pin IO1 changes from low to high state in a Doric DAQ.

% 2019-02-01. Leonardo Molina.
% 2019-08-30. Last modified.
function ttl = loadDoricTTL(filename)
    fid = fopen(filename, 'r');
    parts = textscan(fid, '%f%s%d', 'Delimiter', ',', 'HeaderLines', 1);
    fclose(fid);
    % Get time of transitions from low to high.
    ttlTarget = 'IO1';
    k = ismember(parts{2}, ttlTarget);
    ttl = parts{1}(k);
    state = parts{3}(k);
    rise = diff(state) == +1;
    ttl = ttl([false; rise]);
    % Append stimulation and baseline events.
    ttl = ttl(:)';
end