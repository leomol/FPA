% Read Sabin's annotation file into a dictionary.
% Example:
%   path = 'C:\Users\Molina\Downloads\Trial 1 Annotations.xlsx';
%   events = readAnnotations(path)
%   events('Robot Triggers')
%   events('N PickUpFood')
%   events('entry')
%   events('exit')

% 2023-06-27. Leonardo Molina.
% 2023-06-27. Last modified.
function events = readAnnotations(path)
    data = readcell(path);
    nameLinear = find(cellfun(@ischar, data));
    k = ismember(data(nameLinear), 'Robot Triggers');
    robotTriggerLinear = nameLinear(k);
    [rowStart, colStart] = ind2sub(size(data), robotTriggerLinear);
    [rowEnd, colEnd] = size(data);
    data = data(rowStart:rowEnd, colStart:colEnd);
    header = data(1, 1:2:end);
    data = data(2:end, 1:2:end);
    k = ~cellfun(@isnumeric, data);
    nans = num2cell(NaN(sum(k, 'all'), 1));
    [data{k}] = deal(nans{:});
    data = cell2mat(data);
    events = containers.Map('KeyType', 'char', 'ValueType', 'any');
    for i = 1:numel(header)
        name = header{i};
        column = data(:, i);
        events(name) = column(~isnan(column));
    end
end