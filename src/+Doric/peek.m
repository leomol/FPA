function peek(filename, resamplingFrequency)
    if nargin < 2
        resamplingFrequency = 20;
    end
    datasets = Doric.getDatasets(filename);

    % Group paths by their parent folder.
    groups = containers.Map();
    
    for i = 1:numel(datasets)
        path = datasets{i};
        parts = strsplit(path, '/');
        parentPath = strjoin(parts(1:end-1), '/');
        leafName  = parts{end};
        
        if isKey(groups, parentPath)
            groups(parentPath) = [groups(parentPath), {leafName}];
        else
            groups(parentPath) = {leafName};
        end
    end

    % Process each group to find peaks.
    parents = keys(groups);
    for i = 1:numel(parents)
        parent = parents{i};
        leafNames = groups(parent);
        signalNames = leafNames(~strcmpi(leafNames, 'Time'));
        time = Doric.load(filename, sprintf('%s/Time', parent));
        sourceFrequency = 1 / median(diff(time));
        figure();
        hold('on');
        title(parent);
        for j = 1:numel(signalNames)
            signalName = signalNames{j};
            dataset = sprintf('%s/%s', parent, signalName);
            signal = Doric.load(filename, dataset);
            % Express frequency as a ratio p/q.
            [p, q] = rat(resamplingFrequency / sourceFrequency);
            % Resample: interpolate every p/q/f, upsample by p, filter, downsample by q.
            [signal2, time2] = resample(signal, time, resamplingFrequency, p, q);
            plot(time2, signal2, 'DisplayName', dataset);
        end
        xlabel('time (s)');
        % axis('tight');
        legend('show');
    end
end