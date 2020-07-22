% Loads a mat file exported by LabChart in an intuitive manner.
% 
% [data, comments, units, titles] = loadLabChart(filename)
% 
%   Numerical data and comments are organized in 2D structures where row
%   indices select channels and column indices select blocks:
%     -Numerical data:
%        data(channel, block).time
%        data(channel, block).data
%       
%     -Comments:
%        comments(channel, block).time
%        comments(channel, block).text
%   
%   Units are organized in a cell array where row indices select channels
%   and column indices select blocks:
%      -Data units:
%        units{channel, block}
%     
%   Channel titles are organized in a cell array:
%        titles{channel}
% 
% Example:
%   data = loadLabChart(filename)
%   plot(data(1, 1).time, data(1, 1).data)
% 
% 2020-02-10. Leonardo Molina.
% 2020-07-22. Last modified.
function [data, units, names, comments] = loadLabChart(filename)
    vars = load(filename);
    [nChannels, nBlocks] = size(vars.datastart);
    
    % Sort data.
    data = struct('time', cell(nChannels, nBlocks), 'data', cell(nChannels, nBlocks));
    empty = repmat({zeros(0, 1)}, nChannels, nBlocks);
    [data.time] = deal(empty{:});
    empty = repmat({zeros(0, 1)}, nChannels, nBlocks);
    [data.data] = deal(empty{:});
    offsets = [0; cumsum(seconds(diff(datetime(datevec(vars.blocktimes)))))];
    for c = 1:nChannels
        for b = 1:nBlocks
            frequency = vars.samplerate(c, b);
            if frequency > 0
                selection = vars.datastart(c, b):vars.dataend(c, b);
                n = numel(selection);
                data(c, b).time = offsets(b) + colon(0, n - 1)' / frequency;
                data(c, b).data = vars.data(selection)';
            end
        end
    end
    
    % Sort comments.
    if ~isfield(vars, 'com')
        vars.com = zeros(0, 5);
        vars.comtext = '';
    end
    nComments = size(vars.comtext, 1);
    text = cell(nComments, 1);
    for c = 1:nComments
        text{c} = strtrim(vars.comtext(c, :));
    end
    comments = struct('time', cell(nChannels, nBlocks), 'text', cell(nChannels, nBlocks));
    empty = repmat({zeros(0, 1)}, nChannels, nBlocks);
    [comments.time] = deal(empty{:});
    empty = repmat({{}}, nChannels, nBlocks);
    [comments.text] = deal(empty{:});
    nEntries = size(vars.com, 1);
    for e = 1:nEntries
        channels = vars.com(e, 1);
        b = vars.com(e, 2);
        if channels == -1
            channels = 1:nChannels;
        end
        for c = channels
            timeIndex1 = vars.com(e, 3);
            timeIndex2 = round(timeIndex1 * vars.samplerate(c, b) / vars.tickrate(b));
            time = offsets(b) + timeIndex2 / vars.samplerate(c, b);
            commentIndex = vars.com(e, 5);
            comments(c, b).time = cat(1, comments(c, b).time, time);
            comments(c, b).text = cat(1, comments(c, b).text, text{commentIndex});
        end
    end
    
    % Sort titles.
    names = cell(nChannels, 1);
    for c = 1:nChannels
        names{c} = strtrim(vars.titles(c, :));
    end
    
    % Sort units.
    units = cell(nChannels, nBlocks);
    for c = 1:nChannels
        for b = 1:nBlocks
            if vars.unittextmap(c, b) == -1
                units{c, b} = '';
            else
                units{c, b} = strtrim(vars.unittext(vars.unittextmap(c, b), :));
            end
        end
    end
end