% [labels, start, stop, distance, speed] = CleverSys.load(filename, <sheet>)
% epochs = CleverSys.load(filename, <sheet>)
% Returns a list of event epochs generated by CleverSys.
% 
% Example:
%   filename = 'data\sequential events - frame.xlsx';
%   sheet = 1;
%   [labels, start, finish, distance, speed] = CleverSys.load(filename, sheet);
%   mean1 = mean(speed{1});
%   mean2 = mean(speed{2});
%   plot([mean1, mean2], 'o');
%   xlim([0, 3]);
%   title(sprintf('Compare speeds between "%s" and "%s"', labels{1}, labels{2}));

% 2019-08-20. Leonardo Molina.
% 2023-11-09. Last modified.
function varargout = load(filename, sheet)
    [~, sheetNames, tableFormat] = xlsfinfo(filename);
    
    switch tableFormat
        case {'xlHtml', 'xlCSV'}
            readParameters = {};
        otherwise
            if nargin == 1
                sheet = 1;
            end
            if isnumeric(sheet)
                readParameters = {'Sheet', sheetNames{sheet}};
            else
                readParameters = {'Sheet', sheet};
            end
    end
    
    state = warning('QUERY', 'MATLAB:datetime:AutoConvertStrings');
    warning('OFF', 'MATLAB:datetime:AutoConvertStrings');
    data = readcell(filename, readParameters{:});
    warning(state.state, 'MATLAB:table:ModifiedAndSavedVarnames');
    % Find header.
    [r, c] = find(cellfun(@(c) isequal(c, 'Event'), data));
    titles = data(r, :);
    % Keep data starting from header.
    valid = cellfun(@isstr, data(:, c));
    valid(1:r) = false;
    data = data(valid, :);
    allLabels = data(:, c);
    uniqueLabels = unique(allLabels);
    nUniqueLabels = numel(uniqueLabels);
    % Initialize outputs.
    start = cell(nUniqueLabels, 1);
    finish = cell(nUniqueLabels, 1);
    distance = cell(nUniqueLabels, 1);
    speed = cell(nUniqueLabels, 1);
    epochs = cell(1, 2 * nUniqueLabels);
    % Identify columns.
    startColumn = getColumn(titles, 'From ');
    durationColumn = getColumn(titles, 'Length(');
    distanceColumn = getColumn(titles, 'Dist(');
    speedColumn = getColumn(titles, 'V(');
    % Output is populated according to the columns available.
    if ~isempty(startColumn) && ~isempty(durationColumn)
        m = cell2mat(data(:, [startColumn, durationColumn]));
        m = cumsum(m, 2);
        for u = 1:nUniqueLabels
            label = uniqueLabels{u};
            k = ismember(allLabels, label);
            start{u} = m(k, 1);
            finish{u} = m(k, 2);
            epochs{2 * u - 1} = label;
            epochs{2 * u - 0} = reshape([start{u}, finish{u}]', 1, []);
        end
    end
    if ~isempty(distanceColumn)
        m = cell2mat(data(:, distanceColumn));
        for u = 1:nUniqueLabels
            label = uniqueLabels{u};
            k = ismember(allLabels, label);
            distance{u} = m(k);
        end
    end
    if ~isempty(speedColumn)
        m = cell2mat(data(:, speedColumn));
        for u = 1:nUniqueLabels
            label = uniqueLabels{u};
            k = ismember(allLabels, label);
            speed{u} = m(k);
        end
    end
    if nargout == 1
        varargout = {epochs};
    else
        varargout = {uniqueLabels, start, finish, distance, speed};
    end
end

function column = getColumn(titles, pattern)
    column = find(cellfun(@(title) startsWith(title, pattern), titles));
end