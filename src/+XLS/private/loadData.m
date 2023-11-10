% [data, names] = loadData(filename)
% Load one of several formats of data expected from a recording with either:
% Doric (csv), Multifiber (csv), Inscopix (csv), XLS files, or Axon (abf).
% 
% [data, names] = loadData(filename, sheetName)
% [data, names, sheetName] = loadData(filename, sheetNumber)
% If the file has multiple sheets, use sheetName or sheetNumber to select one.

% 2019-05-07. Leonardo Molina.
% 2021-11-17. Last modified.
function [data, names, sheetName] = loadData(filename, varargin)
    [~, ~, extension] = fileparts(filename);
    if ismember(lower(extension), {'.xls', '.xlsx'})
        [data, names, sheetName] = loadXLS(filename, varargin{:});
    elseif ismember(lower(extension), '.abf')
        [data, ~, names] = ABF.load(filename, varargin{:});
        sheetName = '';
    else
        [data, names] = loadCSV(filename);
        sheetName = '';
    end
end

function [data, names, sheetName] = loadXLS(filename, sheet)
    if nargin == 1
        sheet = 1;
    end
    if isnumeric(sheet)
        [~, sheetNames] = xlsfinfo(filename);
        sheetName = sheetNames{sheet};
    else
        sheetName = sheet;
    end
    state = warning('QUERY', 'MATLAB:table:ModifiedAndSavedVarnames');
    warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames');
    table = readtable(filename, 'Sheet', sheetName);
    warning(state.state, 'MATLAB:table:ModifiedAndSavedVarnames');
    data = table2array(table);
    names = table.Properties.VariableNames;
end

function [data, names] = loadCSV(filename)
    %  <                    , cell1   , cell2 , ... , cellN   >
    %  <Time(s)/Cell Status , accepted,  ...  , ... , accepted>
    %   00.00               , 10.00   ,       , ... , 100.00
    %   00.10               , ..
    % 
    %  < ---                 , channel2 , channel3 , ... , channelN   >
    %   Time(s), channelName1,  channelName2  , ... , channelNameN>
    %   00.00                , 10.00   ,       , ... , 100.00
    %   00.10                , ..
    
    % First line is possibly a header line.
    fid = fopen(filename, 'r');
    tmp = fgetl(fid);
    tmp = textscan(tmp, '%s', 'Delimiter', ',');
    tmp = tmp{1};
	names = tmp;
    if validateNumber(names{1})
        % The file starts with data without a header.
        names = arrayfun(@num2str, -1:numel(names), 'UniformOutput', false);
        nHeaderLines = 0;
    else
        % The file starts with at least one header line.
		tmp = fgetl(fid);
		tmp = textscan(tmp, '%s', 'Delimiter', ',');
		tmp = tmp{1};
        if validateNumber(tmp{1})
            % The file has a single header line.
            nHeaderLines = 1;
        else
            % The file has 2 header lines. Keep the first header.
            nHeaderLines = 2;
        end
    end
    if isempty(names{1})
        % If empty, set name for the first column.
        names{1} = 'time';
    end
    % Next line is data. Detect numeric columns.
    line = fgetl(fid);
    line = textscan(line, '%s', 'Delimiter', ',');
    line = line{1};
    keep = validateNumber(line);
    % Reset carret.
    fseek(fid, 0, 'bof');
    data = textscan(fid, repmat('%f', 1, numel(names)), 'Delimiter', ',', 'HeaderLines', nHeaderLines);
    fclose(fid);
    data = cat(2, data{keep});
    names = names(keep);
    % Remove columns with only NaN values.
    c = all(isnan(data), 1);
    data(:, c) = [];
    names(c) = [];
    % Remove rows with any NaN values.
    r = any(isnan(data), 2);
    data(r, :) = [];
end

function result = validateNumber(numbers)
    result = ~isnan(str2double(numbers)) | contains(numbers, 'nan', 'IgnoreCase', true);
end