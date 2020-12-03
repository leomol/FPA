function [success, messages] = validateEpochRanges(epochRange)
    success = false;
    messages = cell(1, 0);
    if isnumeric(epochRange) && isreal(epochRange) && all(~isnan(epochRange), 'all')
        % Epoch is a vector with numbers that are not complex or nan.
        if mod(numel(epochRange), 2) == 0
            % Epoch range count is even.
            if all(epochRange(2:2:end) >= epochRange(1:2:end), 'all')
                success = true;
            else
                messages{end + 1} = 'epoch ranges must be pairs of increasing values.';
            end
        else
            messages{end + 1} = 'epoch ranges must have an even number of elements.';
        end
    else
        messages{end + 1} = 'epoch range must be a numeric vector.';
    end
end