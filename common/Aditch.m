% LOADADICHT - Loads LabChart's adicht file.
% 
% information = Aditch.getInformation(filename)
% [time, data] = Aditch.getData(filename)

% 2020-02-10. Leonardo Molina.
% 2020-11-20. Last modified.

classdef Aditch
    methods (Static)
        function [time, data] = getData(filename)
            filePointer = getFilePointer(filename);
            nChannels = getNumberOfChannels(filePointer);
            nBlocks = getNumberOfBlocks(filePointer);

            time = repmat({zeros(0, 1)}, nChannels, nBlocks);
            data = repmat({zeros(0, 1)}, nChannels, nBlocks);
            period = NaN(nChannels, nBlocks);
            for c = 1:nChannels
                for b = 1:nBlocks
                    n = getNSamplesInBlock(filePointer, c, b);
                    period(c, b) = getTickPeriod(filePointer, c, b); % !! should offset!?
                    time{c, b} = colon(0, period(c, b), (n - 1) * period(c, b))';
                    data{c, b} = getChannelData(filePointer, c, b, 1, n)';
                end
            end
            closeFileAccessor(filePointer);
        end
        
        function [information, relativeComments, absoluteComments] = getInformation(filename)
            filePointer = getFilePointer(filename);
            nChannels = getNumberOfChannels(filePointer);
            nBlocks = getNumberOfBlocks(filePointer);
            
            names = cell(nChannels, 1);
            units = cell(nChannels, nBlocks);
            samples = NaN(nChannels, nBlocks);
            [dataOffset, blockOffset, period, acquisitionStart] = getBlocksInformation(filePointer);
            for c = 1:nChannels
                names{c} = getChannelName(filePointer, c);
                for b = 1:nBlocks
                    units{c, b} = getUnits(filePointer, c, b);
                    samples(c, b) = getNSamplesInBlock(filePointer, c, b);
                end
            end

            relativeComments = repmat({cell(1, 0)}, nChannels, nBlocks);
            absoluteComments = cell(1, 0);
            for b = 1:nBlocks
                [success, commentPointer] = getCommentPointer(filePointer, b);
                while success
                    [~, tick, comment, channelId] = getComment(commentPointer);
                    if channelId == -1
                        channels = 1:nChannels;
                    else
                        channels = channelId;
                    end
                    for c = channels
                        % Time is relative to the start of the block.
                        offset = dataOffset(c, b) - blockOffset(b);
                        timestamp = tick * period(c, b) + offset;
                        relativeComments{c, b} = cat(2, relativeComments{c, b}, timestamp, comment);
                    end
                    % Time is relative to the start of the first block.
                    % Potential logic error in ADInstruments' SDK:
                    %   If different channels start at different times in a
                    %   given block, and the returned channelId is -1 (aka
                    %   comment is applicable to all channels), the time of
                    %   such comment becomes ambiguous.
                    %   Workaround: Choose channel that starts latest.
                    start = max(dataOffset(channels, b));
                    absoluteComments = cat(2, absoluteComments, timestamp + start, comment);
                    success = advanceComments(commentPointer);
                end
            end

            closeCommentAccessor(commentPointer);
            closeFileAccessor(filePointer);
            
            information.start = acquisitionStart;
            information.names = names;
            information.units = units;
            information.samples = samples;
            information.period = period;
            information.offset = blockOffset;
        end
    end
end

function [dataOffset, blockOffset, period, acquisitionStart] = getBlocksInformation(filePointer)
    % dataOffset:
    %   Time of the first sample in a given channel*block, relative to the
    %   start time of the first block.
    % blockOffset:
    %   Time at which a block was started, relative to the start time of
    %   the first block.
    % 
    % The first data sample in a given channel*block, may occur later than
    % the start of such block.
    
    nChannels = getNumberOfChannels(filePointer);
    nBlocks = getNumberOfBlocks(filePointer);
    dataOffset = zeros(nChannels, nBlocks);
    blockOffset = zeros(1, nBlocks);
    period = NaN(nChannels, nBlocks);
    for c = 1:nChannels
        [dataStartDateB1, acquisitionStart, period(c, 1)] = getBlockInformation(filePointer, c, 1);
        dataOffset(c, 1) = seconds(dataStartDateB1 - acquisitionStart);
        for b = 2:nBlocks
            [dataStartDate, blockStartDate, period(c, b)] = getBlockInformation(filePointer, c, b);
            dataOffset(c, b) = seconds(dataStartDate - dataStartDateB1);
            blockOffset(1, b) = seconds(blockStartDate - acquisitionStart);
        end
    end
end

function [dataStartDate, blockStartDate, period] = getBlockInformation(filePointer, channel, block)
    [~, triggerUnixTime, triggerUnixTimeFraction, triggerMinusStartUnixTicks] = sdk(16, filePointer, zeroIndex(block));
    triggerMinusStartUnixTicks = double(triggerMinusStartUnixTicks);
    blockStartUnix = triggerUnixTime + triggerUnixTimeFraction;
    blockStartDate = datetime(unixToMatlab(blockStartUnix), 'ConvertFrom', 'datenum');
    period = getTickPeriod(filePointer, channel, block);
    dataStartUnix = blockStartUnix - triggerMinusStartUnixTicks * period;
    dataStartDate = datetime(unixToMatlab(dataStartUnix), 'ConvertFrom', 'datenum');
end

function matlabTime = unixToMatlab(unixTime)
    secondsInDay = 86400;
    unixEpoch = 719529;
    matlabTime = unixTime ./ secondsInDay + unixEpoch; 
end
        
function filePointer = getFilePointer(filename)
    [~, filePointer] = sdk(0, [int16(filename), 0]);
    if filePointer == 0
        error('Unable to read file "%s". No such file or directory.', filename);
    end
end

function n = getNumberOfChannels(filePointer)
    [~, n] = sdk(2, filePointer);
    n = double(n);
end

function n  = getNumberOfBlocks(filePointer)
    [~, n] = sdk(1, filePointer);
    n = double(n);
end

function name = getChannelName(filePointer, channel)
    [~, str, n] = sdk(12, filePointer, zeroIndex(channel));
    name = getStringFromOutput(str, n);
end

function units = getUnits(filePointer, channel, block)
    [~, str, n] = sdk(11, filePointer, zeroIndex(block), zeroIndex(channel));
    units = getStringFromOutput(str, n);
end

function n = getNSamplesInBlock(filePointer, channel, block)
    [~, n] = sdk(5, filePointer, zeroIndex(block), zeroIndex(channel));
    n = double(n);
end

function period = getTickPeriod(filePointer, channel, block)
    [~, period] = sdk(4, filePointer, zeroIndex(block), zeroIndex(channel));
end

function data = getChannelData(filePointer, channel, block, startSample, nSamples)
    dataType = int32(0);
    [~, data] = sdk(10, filePointer, zeroIndex(channel), zeroIndex(block), zeroIndex(startSample), int32(nSamples), dataType);
    data = double(data);
end

function closeFileAccessor(filePointer)
    sdk(13, filePointer);
end

function [success, tick, comment, channel] = getComment(commentPointer)
    [result, str, n, tick, channel, ~] = sdk(8, commentPointer);
    if result == 0
        tick = double(tick);
        channel = double(channel);
        comment = getStringFromOutput(str, n);
        success = true;
    else
        comment = '';
        success = false;
    end
end

function success = advanceComments(commentPointer)
    result = sdk(9, commentPointer);
    success = commentParsed(result);
end

function [success, commentPointer] = getCommentPointer(filePointer, block)
    [result, commentPointer] = sdk(6, filePointer, zeroIndex(block));
    success = commentParsed(result);
end

function success = commentParsed(result)
    success = mod(result, 16) ~= 5;
end

function closeCommentAccessor(commentPointer)
    sdk(7, commentPointer);
end

function str = getStringFromOutput(str, n)
    str = char(str(1:n - 1));
end

function output = zeroIndex(input)
    output = int32(input - 1);
end

function varargout = sdk(varargin)
    [varargout{1:nargout}] = resources.sdk_mex(varargin{:});
end