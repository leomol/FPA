% FPA(time, signal, reference, configuration)
%     Plot spontaneous signal from a fiber-photometry experiment based on the
%     values of the configuration. inputFile is a CSV file containing time,
%     signal (465nm) and reference (405nm) columns, or a folder with data
%     recorded with TDT/Synapse.
% 
% Fiber-photometry analysis steps:
%     -Load and resample inputFile.
%     -Low-pass filter an artifact-free portion of the data and fit an exponential
%      decay to correct for photo-bleaching.
%     -Correct for movement artifacts with reference signal.
%     -Compute z-score and low-pass filter to detect peaks.
%     -Compute df/f in a moving time window to normalize traces around peaks.
%     -Compute triggered averages of spontaneous activity grouped by condition/epochs
%      definition.
% 
% configuration is a struct with the following fields (defaults are used for
% missing fields):
%     resamplingFrequency - Resampling frequency (Hz).
%     artifactEpochs - Time epochs (s) to correct via interpolation.
%     zScoreEpochs - Time epochs (s) to include for z-score normalization.
%     bleachingCorrectionEpochs - Time epochs (s) to include for bleaching correction.
%     f0Function - One of @movmean, @movmedian, @movmin.
%     f0Window - Length of the moving window to calculate f0.
%     f1Function - One of @movmean, @movmedian, @movmin, @movstd.
%     f1Window - Length of the moving window to calculate f1.
%     peaksLowpassFrequency - Frequency of lowpass filter to detect peaks.
%     thresholdingFunction - One of @mad, @std.
%     thresholdFactor - Thresholding cut-off.
%     conditionEpochs - Epochs that involve different conditions: {'epoch1', [start1, end1, start2, end2, ...], 'epoch2', ...}
%     triggeredWindow - Length of the time window around each peak.
% 
% See source code for default values.
% Units for time and frequency are seconds and hertz respectively.
% 
% Normalization:
%     df/f is calculated as (f - f0) / f1, where f0 and f1 are computed for
%     each time point using function over a moving window with the given size (s).
%     Recipe for df/f:
%         configuration.f0Function = @movmean;
%         configuration.f1Function = @movmean;
%     Recipe for z-score:
%         configuration.f0Function = @movmean;
%         configuration.f1Function = @movstd;
% 
% See FPAexamples

% 2019-02-01. Leonardo Molina.
% 2019-08-30. Last modified.
function results = FPA(time, signal, reference, configuration)
    addpath(fullfile(fileparts(mfilename('fullpath')), 'common'));
    if nargin == 1
        configuration = struct();
    end
    
    configuration = setDefault(configuration, 'resamplingFrequency', 50);
    configuration = setDefault(configuration, 'artifactEpochs', []);
    configuration = setDefault(configuration, 'zScoreEpochs', [-Inf, Inf]);
    configuration = setDefault(configuration, 'bleachingCorrectionEpochs', [-Inf, Inf]);
    configuration = setDefault(configuration, 'f0Function', @movmean);
    configuration = setDefault(configuration, 'f0Window', 10);
    configuration = setDefault(configuration, 'f1Function', @movstd);
    configuration = setDefault(configuration, 'f1Window', 10);
    configuration = setDefault(configuration, 'peaksLowpassFrequency', 0.2);
    configuration = setDefault(configuration, 'thresholdingFunction', @std);
    configuration = setDefault(configuration, 'thresholdFactor', 0.10);
    configuration = setDefault(configuration, 'conditionEpochs', {'Condition A', [-Inf, Inf]});
    configuration = setDefault(configuration, 'triggeredWindow', 10);
    sourceFrequency = 1 / mean(diff(time));
    
    % Resample to target frequency.
    [p, q] = rat(configuration.resamplingFrequency / sourceFrequency);
    reference = resample(reference, time, configuration.resamplingFrequency, p, q);
    signal = resample(signal, time, configuration.resamplingFrequency, p, q);
    time = linspace(time(1), numel(signal) / configuration.resamplingFrequency, numel(signal))';
    
    % Remove artifacts by interpolating data.
    id = transpose(1:numel(time));
    artifactId = setdiff(id, time2id(time, configuration.artifactEpochs));
    id2 = id(artifactId);
    r2 = reference(artifactId);
    s2 = signal(artifactId);
    r2 = interp1(id2, r2, id);
    s2 = interp1(id2, s2, id);
    
    % Indexed epoch range.
    bleachingCorrectionId = time2id(time, configuration.bleachingCorrectionEpochs);
    
    % Remove peaks for bleaching correction.
    filterOrder = 12;
    bleachingLowpassFrequency = 0.1;
    bleachingFilter = designfilt('lowpassiir', 'HalfPowerFrequency', bleachingLowpassFrequency, 'SampleRate', configuration.resamplingFrequency, 'DesignMethod', 'butter', 'FilterOrder', filterOrder);
    
    % Fitting function.
    rLowpass = filtfilt(bleachingFilter, r2);
    sLowpass = filtfilt(bleachingFilter, s2);
    rFit = fit(time(bleachingCorrectionId), rLowpass(bleachingCorrectionId), fittype('exp1'));
    sFit = fit(time(bleachingCorrectionId), sLowpass(bleachingCorrectionId), fittype('exp1'));
    rBleaching = rFit(time);
    sBleaching = sFit(time);
    rPost = r2 ./ rBleaching;
    sPost = s2 ./ sBleaching;
    
    % Normalize to reference signal.
    f = sPost ./ rPost;
    
    % Standardize.
    zScoreId = time2id(time, configuration.zScoreEpochs);
    z = (f - mean(f(zScoreId))) / std(f(zScoreId));
    
    % Normalization factor.
    f0 = configuration.f0Function(f, nanmin(round(configuration.f0Window * configuration.resamplingFrequency), numel(f)));
    f1 = configuration.f1Function(f, nanmin(round(configuration.f1Window * configuration.resamplingFrequency), numel(f)));
    df = f - f0;
    dff = df ./ f1;
    
    % Find peaks.
    filterOrder = 12;
    peaksFilter = designfilt('lowpassiir', 'HalfPowerFrequency', configuration.peaksLowpassFrequency, 'SampleRate', configuration.resamplingFrequency, 'DesignMethod', 'butter', 'FilterOrder', filterOrder);
    signalLowPass = filtfilt(peaksFilter, z);
    threshold = mean(signalLowPass) + configuration.thresholdFactor * configuration.thresholdingFunction(signalLowPass);
    [~, peaksId] = findpeaks(signalLowPass, 'MinPeakHeight', threshold);
    [~, valleysId] = findpeaks(-signalLowPass);
    
    % Find the actual peak around a time window compatible with the filter.
    w = ceil(0.5 * configuration.resamplingFrequency / configuration.peaksLowpassFrequency);
    if mod(w, 2) ~= 0
        w = w + 1;
    end
    nSamples = numel(time);
    for i = 1:numel(peaksId)
        range = max(1, peaksId(i) - w / 2):min(nSamples, peaksId(i) + w / 2);
        [~, k] = max(+z(range));
        peaksId(i) = k + range(1) - 1;
        [~, k] = max(-z(range));
        valleysId(i) = k + range(1) - 1;
    end
    
    % Keep highest peak within a period.
    keep = true(size(peaksId));
    for i = 1:numel(valleysId) - 1
        k = peaksId > valleysId(i) & peaksId < valleysId(i);
        [~, largest] = max(z(k));
        m = false(sum(k), 1);
        m(largest) = true;
        keep(k) = m;
    end
    peaksId = peaksId(keep);
    
    % Split peaks/traces by conditions.
    allIds = colon(1, nSamples)';
    conditionBool = cell(size(configuration.conditionEpochs));
    % Number of samples in a window (force odd number).
    window = round(configuration.triggeredWindow * configuration.resamplingFrequency);
    if mod(window, 2) == 0
        window = window - 1;
    end
    % Index template to apply on each peak.
    windowTemplate = -(window - 1) / 2:(window - 1) / 2;
    % Filter out out of range traces.
    peaksId2 = peaksId(peaksId > (window - 1) / 2 & peaksId + (window - 1) / 2 < nSamples);
    % Index of all triggered traces, one per row.
    if numel(peaksId2) > 0
        triggeredId = bsxfun(@plus, windowTemplate, peaksId2);
    else
        triggeredId = [];
    end
    % Label each row according to the conditions.
    conditionIds = zeros(1, numel(peaksId2));
    nEpochs = numel(configuration.conditionEpochs) / 2;
    for id = 1:2:2 * nEpochs
        conditionBool(id:id + 1) = [configuration.conditionEpochs(id), ismember(allIds, time2id(time, configuration.conditionEpochs{id + 1}))];
        k = conditionBool{id + 1}(peaksId2);
        conditionIds(k) = (id + 1) / 2;
    end
    
    % Plots.
    cmap = lines();
    
    % Plot raw signal and bleaching.
    figure();
    ax.pre = axes('XTick', []);
    hold(ax.pre, 'all');
    ax.pre.Position = [ax.pre.Position(1), 0.7, ax.pre.Position(3), 0.3];
    plot(ax.pre, time, signal);
    plot(ax.pre, time, sBleaching, '--');
    legend(ax.pre, {'Signal', 'Bleaching'});
    
    % Plot corrected signal, low-pass filtered signal and peaks.
    ax.post = axes('XTick', []);
    hold(ax.post, 'all');
    ax.post.Position = [ax.post.Position(1), 0.40, ax.post.Position(3), 0.25];
    plot(ax.post, time, z);
    plot(ax.post, time, signalLowPass, 'k');
    plot(ax.post, time([1, end]), threshold([1, 1]), 'k--');
    plot(ax.post, time(peaksId), z(peaksId), 'ro');
    legend(ax.post, {'z-score', 'low-pass', 'threshold', 'peaks'});
    
    % Plot traces around peaks.
    ax.dff = axes();
    xlabel(ax.dff, 'Time (s)');
    hold(ax.dff, 'all');
    ax.dff.Position = [ax.dff.Position(1), 0.1, ax.dff.Position(3), 0.25];
    t2 = [time(triggeredId) NaN(size(triggeredId, 1), 1)]';
    y2 = [dff(triggeredId) NaN(size(triggeredId, 1), 1)]';
    t3 = time(peaksId);
    y3 = dff(peaksId);
    ylims = [min(y2(:)), max(y2(:))];
    for e = 1:nEpochs
        epochName = configuration.conditionEpochs{2 * e - 1};
        [faces, vertices] = patchEpochs(configuration.conditionEpochs{2 * e}, ylims(1), ylims(2));
        patch('Faces', faces, 'Vertices', vertices, 'FaceColor', cmap(e, :), 'EdgeColor', 'none', 'FaceAlpha', 0.50, 'DisplayName', sprintf('%s', epochName));
    end
    plot(ax.dff, t2(:), y2(:), 'DisplayName', 'df/f');
    plot(ax.dff, t3, y3, 'Color', 'r', 'LineStyle', 'none', 'Marker', 'o', 'DisplayName', 'peaks');
    legend(ax.dff, 'show');
    
    linkaxes([ax.pre, ax.post, ax.dff], 'x');
    
    % Plot triggered average.
    figure('name', 'FPA');
    ax.trigger = axes();
    hold(ax.trigger, 'all');
    uIds = unique(conditionIds);
    uIds = uIds(uIds > 0);
    triggerT = windowTemplate / configuration.resamplingFrequency;
    nPeaks = zeros(size(uIds));
    semAtZero = zeros(size(uIds));
    stdPerCondition = zeros(size(uIds));
    for i = 1:numel(uIds)
        id = uIds(i);
        nPeaks(i) = sum(conditionIds == id);
        triggeredData = dff(triggeredId(conditionIds == id, :));
        triggeredData = reshape(triggeredData, numel(triggeredData) / window, window);
        triggerMean = mean(triggeredData, 1);
        triggerSem = std(triggeredData, [], 1) / sqrt(size(triggeredData, 1));
        semAtZero(i) = triggerSem(ceil(window / 2));
        stdPerCondition(i) = std(z(conditionBool{2 * i}));
        h1 = plot(triggerT, triggerMean);
        vertices = [triggerT; triggerMean + triggerSem / 2];
        vertices = cat(2, vertices, [fliplr(triggerT); fliplr(triggerMean - triggerSem / 2)])';
        faces = 1:2 * window;
        patch('Faces', faces, 'Vertices', vertices, 'FaceColor', h1.Color, 'EdgeColor', 'none', 'FaceAlpha', 0.10);
    end
    legend(gca, [configuration.conditionEpochs(2 * uIds - 1); arrayfun(@(i) ['n=' sprintf('%i', nPeaks(i)) ', SEM=' sprintf('%.4f', semAtZero(i))], 1:numel(uIds), 'UniformOutput', false)]);
    h = title('Triggered average for each condition \pm SEM');
    set(h, 'Interpreter', 'tex');
    xlabel('Time (s)');
    ylabel('df/f');
    
    % Plot stats on traces.
    figure('name', 'FPA');
    counts = zeros(nEpochs, 1);
    epochIds = zeros(1, 0);
    for e = 1:nEpochs
        k = time2id(time, configuration.conditionEpochs{2 * e});
        epochIds = cat(2, epochIds, k);
        counts(e) = numel(k);
    end
    labelIds = zeros(1, sum(counts));
    labelIds(cumsum(counts(1:end-1)) + 1) = 1;
    labelIds = cumsum(labelIds) + 1;
    keep = ismember(labelIds, uIds);
    labelIds = labelIds(keep);
    epochIds = epochIds(keep);
    labels = arrayfun(@(i) sprintf('"%s" STD:%.2f', configuration.conditionEpochs{2 * i - 1}, stdPerCondition(i)), 1:numel(uIds), 'UniformOutput', false);
    boxplot(z(epochIds), labelIds, 'Labels', labels);
    ylabel('z-score');
    title('Stats on z-scored traces for each epoch');
    
    axis([ax.pre, ax.dff, ax.trigger], 'tight');
    
    results.peaksId = peaksId;
    results.dff = dff;
end

function configuration = setDefault(configuration, fieldname, value)
    if ~isfield(configuration, fieldname)
        configuration.(fieldname) = value;
    end
end