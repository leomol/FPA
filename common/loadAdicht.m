% LOADADICHT - Loads LabChart's adicht file.
% 
% [time, data, units, names] = loadAdicht(filename)
% 
%   Data are organized in 2D cells where row indices select channels and
%   column indices select blocks:
%        time{channel, block}
%        data{channel, block}
%   
%   Units are organized in a cell array where row indices select channels
%   and column indices select blocks:
%        units{channel, block}
%     
%   Channel names are organized in a cell array:
%        titles{channel}
% 
% Example:
%   data = loadLabChart(filename)
%   plot(time{1, 1}, data{1, 1})
% 
% See LOADLABCHART.
% 
% 2020-02-10. Leonardo Molina.
% 2020-10-22. Last modified.
function [time, data, units, names] = loadAdicht(filename)
    pointer = openFile(filename);
    
    nChannels = getNumberOfChannels(pointer);
    nBlocks = getNumberOfRecords(pointer);
    
    time = repmat({zeros(0, 1)}, nChannels, nBlocks);
    data = repmat({zeros(0, 1)}, nChannels, nBlocks);
    units = cell(nChannels, nBlocks);
    names = cell(nChannels, 1);
    for c = 1:nChannels
        names{c} = getChannelName(pointer, c);
        for b = 1:nBlocks
            units{c, b} = getUnits(pointer, b, c);
            n = getNSamplesInRecord(pointer, b, c);
            dt = getTickPeriod(pointer, b, c);
            time{c, b} = colon(0, dt, (n - 1) * dt)';
            data{c, b} = getChannelData(pointer, b, c, 1, n)';
        end
    end
    
    closeFile(pointer);
end

function pointer = openFile(filename)
    [~, pointer] = sdk(0, [int16(filename), 0]);
    if pointer == 0
        error('Unable to read file "%s". No such file or directory.', filename);
    end
end

function n = getNumberOfChannels(pointer)
    [~, n] = sdk(2, pointer);
    n = double(n);
end

function n  = getNumberOfRecords(pointer)
    [~, n] = sdk(1, pointer);
    n = double(n);
end

function name = getChannelName(pointer, channel)
    [~, str, n] = sdk(12, pointer, c0(channel));
    name = getStringFromOutput(str, n);
end

function units = getUnits(pointer, record, channel)
    [~, str, n] = sdk(11, pointer, c0(record), c0(channel));
    units = getStringFromOutput(str, n);
end

function n = getNSamplesInRecord(pointer, record, channel)
    [~, n] = sdk(5, pointer, c0(record), c0(channel));
    n = double(n);
end

function period = getTickPeriod(pointer, record, channel)
    [~, period] = sdk(4, pointer, c0(record), c0(channel));
end

function data = getChannelData(pointer, record, channel, startSample, nSamples)
    dataType = c(0);
    [~, data] = sdk(10, pointer, c0(channel), c0(record), c0(startSample), c(nSamples), dataType);
    data = double(data);
end

function closeFile(pointer)
    sdk(13, pointer);
end

function str = getStringFromOutput(str, n)
    str = char(str(1:n - 1));
end

function output = c0(input)
    output = int32(input - 1);
end

function output = c(input)
    output = int32(input);
end

function varargout = sdk(varargin)
    [varargout{1:nargout}] = resources.sdk_mex(varargin{:});
end