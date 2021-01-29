% [data, names] = loadData(filename)
% Load one of several formats of data expected from a recording with either:
% Doric (csv), Multifiber (csv), Inscopix (csv), XLS files, or Axon (abf).
% 
% [data, names] = loadData(filename, sheetName)
% [data, names, sheetName] = loadData(filename, sheetNumber)
% If the file has multiple sheets, use sheetName or sheetNumber to select one.

% 2019-05-07. Leonardo Molina.
% 2021-01-28. Last modified.
function [data, names, sheetName] = loadData(filename, varargin)
    [~, ~, extension] = fileparts(filename);
    if ismember(lower(extension), {'.xls', '.xlsx'})
        [data, names, sheetName] = loadXLS(filename, varargin{:});
    elseif ismember(lower(extension), '.abf')
        [data, ~, names] = loadABF(filename, varargin{:});
        sheetName = '';
    else
        [data, names] = loadCSV(filename, varargin{:});
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

function [data, names] = loadCSV(filename, nRows)
    %  <                    , cell1   , cell2 , ... , cellN   >
    %  <Time(s)/Cell Status , accepted,  ...  , ... , accepted>
    %   00.00               , 10.00   ,       , ... , 100.00
    %   00.10               , ..
    % 
    %  < ---                 , channel2 , channel3 , ... , channelN   >
    %   Time(s), channelName1,  channelName2  , ... , channelNameN>
    %   00.00                , 10.00   ,       , ... , 100.00
    %   00.10                , ..
    
    if nargin == 1
        nRows = Inf;
    end
    
    fid = fopen(filename, 'r');
    tmp = fgetl(fid);
    tmp = textscan(tmp, '%s', 'Delimiter', ',');
    tmp = tmp{1};
	names = tmp;
    if isnan(str2double(names{1}))
		tmp = fgetl(fid);
		tmp = textscan(tmp, '%s', 'Delimiter', ',');
		tmp = tmp{1};
		if isnan(str2double(tmp{1}))
			nHeaderLines = 2;
			names = tmp;
		else
			nHeaderLines = 1;
		end
    else
        names = arrayfun(@num2str, -1:numel(names), 'UniformOutput', false);
        nHeaderLines = 0;
    end
    if isempty(names{1})
        names{1} = 'time';
    end
    keep = fgetl(fid);
    keep = textscan(keep, '%s', 'Delimiter', ',');
    keep = keep{1};
    fseek(fid, 0, 'bof');
    if ismember({'accepted'}, keep)
        nHeaderLines = nHeaderLines + 1;
        keep = [true; ismember(keep(2:end), 'accepted')];
        data = textscan(fid, repmat('%f', 1, numel(names)), 'Delimiter', ',', 'HeaderLines', nHeaderLines);
    else
        keep = ~isnan(str2double(keep));
        format = cell(1, numel(keep));
        [format{ keep}] = deal('%f');
        [format{~keep}] = deal('%s');
        format = cat(2, format{:});
        data = textscan(fid, format, nRows, 'Delimiter', ',', 'HeaderLines', nHeaderLines);
    end
    data = cat(2, data{keep});
    names = names(keep);
    fclose(fid);
end