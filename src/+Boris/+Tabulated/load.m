% [labels, start, stop] = Boris.Tabulated.load(filename)
% Returns a list of events encoded with BORIS.
% Output:
%   labels: 'behavior1', 'behavior2', ...
%   start: [start1, start2, ...];
%   stop: [stop1, stop2, ...]
% 
% epochs = Boris.Tabulated.load(filename)
% Output:
%   {'behavior1', [start1, stop1, start2, stop2, ...], 'behavior2', [start1, stop1, start2, stop2, ...]}

% 2021-02-26. Leonardo Molina.
% 2024-01-19. Last modified.
function varargout = load(filename)
    % Data is separated by commas or tabs.
    fid = fopen(filename, 'r');
    line = fgetl(fid);
    delimiter = ',';
    header = strsplit(line, delimiter);
    if numel(header) == 1
        delimiter = '\t';
        header = strsplit(line, delimiter);
    end
    % Read all columns as text.
    format = repmat('%s', size(header));
    data = textscan(fid, format, 'Delimiter', delimiter);
    fclose(fid);
    % Parse target columns.
    labelsColumn = ismember(header, 'Behavior');
    startColumn = ismember(header, 'Start (s)');
    stopColumn = ismember(header, 'Stop (s)');
    labels = data{labelsColumn};
    start = str2double(data{startColumn});
    stop = str2double(data{stopColumn});

    % Number of outputs determines output format.
    if nargout == 1
        uLabels = unique(labels);
        nLabels = numel(uLabels);
        epochs = cell(1, 2 * nLabels);
        epochs(1:2:end) = uLabels;
        for u = 1:nLabels
            label = uLabels{u};
            k = ismember(labels, label);
            epochs{2 * u} = [start(k), stop(k)]';
        end
        varargout = {epochs};
    else
        varargout = {labels, start, stop};
    end
end