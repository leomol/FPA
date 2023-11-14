% DLC.load
% Read DLC in csv format into a MATLAB table.

% 2022-07-13. Leonardo Molina.
% 2023-11-09. Last modified.
function data = load(path)
    fid = fopen(path, 'rt');
    % Ignore first line.
    fgetl(fid);
    % Get header.
    line = fgetl(fid);
    columns = strsplit(line, ',');
    nColumns = numel(columns);
    % Reset carret.
    fseek(fid, 0, 'bof');
    % Read all, drop frame index.
    data = textscan(fid, repmat('%f', 1, nColumns), 'Delimiter', ',', 'HeaderLines', 3);
    fclose(fid);
    data = [data{2:end}];
    % Append x, y, p to header.
    header = cellfun(@(part) cellfun(@(suffix) sprintf('%s_%s', part, suffix), {'x', 'y', 'p'}, 'UniformOutput', false), columns(2:3:end), 'UniformOutput', false);
    header = [header{:}];
    % Return table.
    data = array2table(data, 'VariableNames', header);
end