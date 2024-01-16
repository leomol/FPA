% [labels, start, stop] = Boris.Aggregated.load(filename)
% Returns a list of events encoded with BORIS.
% Output:
%   labels: 'behavior1', 'behavior2', ...
%   start: [start1, start2, ...];
%   stop: [stop1, stop2, ...]
% 
% epochs = Boris.Aggregated.load(filename)
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
    format = repmat({'%s'}, size(header));
    % Some columns are expected to be numeric.
    [format{ismember(header, {'Time', 'Total length', 'FPS'})}] = deal('%f');
    format = [format{:}];
    % Target columns.
    [~, columns] = intersect(header, {'Behavior', 'Behavior type', 'Time'});
    % Reset cursor.
    fseek(fid, 0, 'bof');
    % Read, sort columns, and assign.
    data = textscan(fid, format, 'Delimiter', delimiter, 'HeaderLines', 1);
    data = data(columns);
    time = cat(1, data{3});
    labels = data{1};
    status = (data{2} == "STOP") + 1;
    
    % Sort data so that every behavior starts and stops in consecutive rows.
    fseek(fid, 0, 'bof');
    format = repmat('%s', size(header));
    % Read time as text to find the largest decimal count.
    data = textscan(fid, format, 'Delimiter', delimiter, 'HeaderLines', 1);
    timeText = data{columns(3)};
    decimalsText = regexp(timeText, '\.(\d+)', 'tokens', 'once');
    decimalsText = [decimalsText{:}];
    decimalsCount = cellfun(@numel, decimalsText);
    % Define function for padding with zeros.
    padLeft = @(s, n) sprintf('%0*s', n, s);
    fixRight = @(x, n) sprintf('%.*f', n, x);
    pad = @(x, l, r) padLeft(fixRight(x, r), l + r + 1);
    % Pad left according to largest integer.
    nLeft = numel(num2str(ceil(max(time))));
    % Pad right accorting to largest decimal count.
    nRight = max(decimalsCount);
    % Create a time string that can be sorted alphabetically.
    timeUID = arrayfun(@(x) pad(x, nLeft, nRight), time, 'UniformOutput', false);
    % Create a row string that can be sorted alphabetically by label > time > status.
    uid = [labels, timeUID, num2cell(num2str(status))]';
    uid = strcat(uid(1, :), uid(2, :), uid(3, :));
    [~, order] = sort(uid);
    start = time(order(1:2:end));
    stop = time(order(2:2:end));
    stop(end + 1:numel(start)) = Inf;
    labels = labels(order(1:2:end));
    
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
    fclose(fid);
end