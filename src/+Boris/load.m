% [labels, start, stop] = Boris.load(filename)
% Returns a list of events encoded with BORIS.
% Output:
%   labels: 'behavior1', 'behavior2', ...
%   start: [start1, start2, ...];
%   stop: [stop1, stop2, ...]
% 
% epochs = Boris.load(filename)
% Output:
%   {'behavior1', [start1, stop1, start2, stop2, ...], 'behavior2', [start1, stop1, start2, stop2, ...]}

% 2021-02-26. Leonardo Molina.
% 2024-01-19. Last modified.
function varargout = load(filename)
    % Data is separated by commas.
    fid = fopen(filename, 'r');
    line = fgetl(fid);
    fclose(fid);
    header = strsplit(line, '[,\t]', 'DelimiterType', 'RegularExpression');
    varargout = cell(1, nargout);
    if ismember('Start (s)', header)
        [varargout{:}] = deal(Boris.Tabulated.load(filename));
    else
        [varargout{:}] = deal(Boris.Aggregated.load(filename));
    end
end