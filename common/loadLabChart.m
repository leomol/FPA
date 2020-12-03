% LOADLABCHART - Loads a dot-mat file exported from LabChart.
% 
% [time, data, units, comments, names] = loadLabChart(filename)
% 
%   Most data are organized in 2D cells where row indices select channels
%   and column indices select blocks:
%        time{channel, block}
%        data{channel, block}
%        units{channel, block}
%        comments{channel, block}
%     
%   Channel names are organized in a cell array:
%        names{channel}

% 2020-02-10. Leonardo Molina.
% 2020-11-19. Last modified.
function [time, data, units, comments, names] = loadLabChart(filename)
    vars = load(filename);
    [nChannels, nBlocks] = size(vars.datastart);
    
    % Sort data.
    time = repmat({zeros(0, 1)}, nChannels, nBlocks);
    data = repmat({zeros(0, 1)}, nChannels, nBlocks);
    offsets = [0; cumsum(seconds(diff(datetime(datevec(vars.blocktimes)))))];
    for c = 1:nChannels
        for b = 1:nBlocks
            frequency = vars.samplerate(c, b);
            if frequency > 0
                selection = vars.datastart(c, b):vars.dataend(c, b);
                n = numel(selection);
                time{c, b} = offsets(b) + colon(0, n - 1)' / frequency;
                data{c, b} = vars.data(selection)';
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
    
    comments = repmat({cell(1, 0)}, nChannels, nBlocks);
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
            timestamp = offsets(b) + timeIndex2 / vars.samplerate(c, b);
            commentIndex = vars.com(e, 5);
            comments{c, b} = cat(2, comments{c, b}, timestamp, text{commentIndex});
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