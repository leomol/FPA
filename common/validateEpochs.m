function [success, messages] = validateEpochs(epochs)
    success = false;
    messages = cell(1, 0);
    if iscell(epochs)
        % Content is cell.
        if mod(numel(epochs), 2) == 0
            % Cell count is even.
            epochNames = epochs(1:2:end);
            if all(cellfun(@isstr, epochNames) | cellfun(@isstring, epochNames))
                % Epoch names are chars or strings.
                uniqueEpochNames = unique(epochNames);
                if numel(uniqueEpochNames) == numel(epochNames)
                    % Epoch names are unique.
                    epochRanges = epochs(2:2:end);
                    nEpochs = numel(epochRanges);
                    success = true;
                    for i = 1:nEpochs
                        [s, m] = validateEpochRanges(epochRanges{i});
                        if s == false
                            % Epoch range is invalid.
                            messages = cat(2, messages, m);
                            success = false;
                            break;
                        end
                    end
                else
                    messages{end + 1} = 'epoch names must be unique.';
                end
            else
                messages{end + 1} = 'epoch names must be chars or strings.';
            end
        else
            messages{end + 1} = 'epoch definitions must have an even number of elements.';
        end
    else
        messages{end + 1} = 'epoch definitions must be a cell array.';
    end
end