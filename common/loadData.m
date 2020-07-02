% [data, names] = loadData(filename)
% Load one of several formats of data expected from a recording with either:
% Doric (csv), Inscopix (xls), or Axon (abf).
% 
% [data, names] = loadData(filename, sheetName)
% [data, names, sheetName] = loadData(filename, sheetNumber)
% If the file has multiple sheets, use sheetName or sheetNumber to select one.

% 2019-05-07. Leonardo Molina.
% 2020-06-08. Last modified.
function [data, names, sheetName] = loadData(filename, sheet)
    [~, ~, extension] = fileparts(filename);
    if ismember(lower(extension), {'.xls', '.xlsx'})
        if nargin == 1
            sheet = 1;
        end
        [data, names, sheetName] = loadXLS(filename, sheet);
    elseif ismember(lower(extension), '.abf')
        [data, ~, names] = loadABF(filename);
    else
        [data, names] = loadCSV(filename);
    end
end

function [data, names, sheetName] = loadXLS(filename, sheet)
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
    names{1} = 'time';
    keep = fgetl(fid);
    keep = textscan(keep, '%s', 'Delimiter', ',');
    keep = keep{1};
    fseek(fid, 0, 'bof');
    if ismember({'accepted'}, keep)
        nHeaderLines = nHeaderLines + 1;
        keep = [true; ismember(keep(2:end), 'accepted')];
        data = textscan(fid, repmat('%f', 1, numel(names)), 'Delimiter', ',', 'HeaderLines', nHeaderLines);
    else
        keep = true(size(names));
        data = textscan(fid, repmat('%f', 1, numel(names)), 'Delimiter', ',', 'HeaderLines', nHeaderLines);
    end
    data = cat(2, data{:});
    data = data(:, keep);
    names = names(keep);
    fclose(fid);
end