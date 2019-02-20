% FPA(inputFile, configuration)
%     Plot spontaneous signal from a fiber-photometry experiment based on the
%     values of the configuration. inputFile is a CSV file containing time,
%     signal and reference columns.
% 
% Fiber photometry analysis:
%     -Load and resample inputFile.
%     -Low-pass filter an artifact-free portion of the data and fit an exponential
%      decay to correct for bleaching.
%     -Correct for movement artifacts with reference signal.
%     -Compute z-score over the whole recording and low-pass filter to detect peaks.
%     -Compute df/f in a moving time window to normalize traces around peaks.
%     -Compute triggered averages of spontaneous activity grouped by condition/epochs
%      definition.
% 
% configuration is a struct with the following fields (defaults are used for
% missing fields):
%     timeTitle - Exact title name of time column.
%     signalTitle - Exact title name of time column.
%     referenceTitle - Exact title name of time column.
%     resamplingFrequency - Resampling frequency (Hz).
%     bleachingCorrectionEpochs - Time epochs (s) to include for bleaching correction.
%     f0Function - One of @movmean, @movmedian, @movmin.
%     f0Window - Length of the moving window to calculate f0.
%     f1Function - One of @movmean, @movmedian, @movmin, @movstd.
%     f1Window - Length of the moving window to calculate f1.
%     peaksLowpassFrequency - Frequency of lowpass filter to detect peaks.
%     thresholdingFunction - One of @mad, @std.
%     thresholdFactor - Thresholding cut-off.
%     conditionEpochs - Epochs that involving different conditions: {'epoch1', [start1, end1, start2, end2, ...], 'epoch2', ...}
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
% Example 1:
%     inputFile = 'data/tf_CP_m307_d2_2.csv';
%     FPA(inputFile);
% 
% Example 2:
%     inputFile = 'data/tf_CP_m307_d2_2.csv';
%     configuration.timeTitle = 'Time(s)';
%     configuration.signalTitle = 'AIn-2 - Demodulated(Lock-In)';
%     configuration.referenceTitle = 'AIn-1 - Demodulated(Lock-In)';
%     configuration.resamplingFrequency = 50;
%     configuration.bleachingCorrectionEpochs = [-Inf, 600, 960, Inf];
%     configuration.conditionEpochs = {'Condition A', [-Inf, 600], 'Condition B', [650, Inf]};
%     configuration.triggeredWindow = 10;
%     configuration.f0Function = @movmean;
%     configuration.f0Window = 10;
%     configuration.f1Function = @movstd;
%     configuration.f1Window = 10;
%     configuration.peaksLowpassFrequency = 0.2;
%     configuration.thresholdingFunction = @std;
%     configuration.thresholdFactor = 0.10;
%     FPA(inputFile, configuration);

% 2019-02-01. Leonardo Molina.
% 2019-02-20. Last modified.
function FPA(inputFile, configuration)
    if nargin == 1
        configuration = struct();
    end
    configuration = setDefault(configuration, 'timeTitle', 'Time(s)');
    configuration = setDefault(configuration, 'signalTitle', 'AIn-2 - Demodulated(Lock-In)');
    configuration = setDefault(configuration, 'referenceTitle', 'AIn-1 - Demodulated(Lock-In)');
    configuration = setDefault(configuration, 'resamplingFrequency', 50);
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
    
    % Processing.
    % Load data from file.
    fid = fopen(inputFile, 'r');
    header = fgetl(fid);
    header = strsplit(header, ',');
    header = header(~cellfun(@isempty, header));
    fclose(fid);
    data = csvread(inputFile, 1, 0);
    t = data(:, find(ismember(header, configuration.timeTitle), 1));
    r = data(:, find(ismember(header, configuration.signalTitle), 1));
    s = data(:, find(ismember(header, configuration.referenceTitle), 1));
    sourceFrequency = 1 / mean(diff(t));
    
    % Resample to target frequency.
    [p, q] = rat(configuration.resamplingFrequency / sourceFrequency);
    r = resample(r, t, configuration.resamplingFrequency, p, q);
    s = resample(s, t, configuration.resamplingFrequency, p, q);
    t = linspace(t(1), numel(s) / configuration.resamplingFrequency, numel(s))';
    
    % Indexed epoch range.
    bleachingCorrectionId = time2id(t, configuration.bleachingCorrectionEpochs);
    
    % Remove peaks for bleaching correction.
    filterOrder = 12;
    bleachingLowpassFrequency = 0.1;
    bleachingFilter = designfilt('lowpassiir', 'HalfPowerFrequency', bleachingLowpassFrequency, 'SampleRate', configuration.resamplingFrequency, 'DesignMethod', 'butter', 'FilterOrder', filterOrder);
    
    % Fitting function.
    rLowpass = filtfilt(bleachingFilter, r);
    sLowpass = filtfilt(bleachingFilter, s);
    rFit = fit(t(bleachingCorrectionId), rLowpass(bleachingCorrectionId), fittype('exp1'));
    sFit = fit(t(bleachingCorrectionId), sLowpass(bleachingCorrectionId), fittype('exp1'));
    rBleaching = rFit(t);
    sBleaching = sFit(t);
    rPost = r ./ rBleaching;
    sPost = s ./ sBleaching;
    
    % Normalize to reference signal.
    f = sPost ./ rPost;
    % Standardize.
    z = (f - mean(f)) / std(f);
    
    % Normalization factor.
    f0 = configuration.f0Function(f, min(round(configuration.f0Window * configuration.resamplingFrequency), numel(f)));
    f1 = configuration.f1Function(f, min(round(configuration.f1Window * configuration.resamplingFrequency), numel(f)));
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
    nSamples = numel(t);
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
    windowTemplateIds = -(window - 1) / 2:(window - 1) / 2;
    % Filter out traces that would have gone out of bounds.
    peaksId2 = peaksId(peaksId > (window - 1) / 2 & peaksId + (window - 1) / 2 < nSamples);
    % Index of all triggered traces, one per row.
    if numel(peaksId2) > 0
        triggeredId = bsxfun(@plus, windowTemplateIds, peaksId2);
    else
        triggeredId = [];
    end
    % Label each row according to the conditions.
    conditionIds = zeros(1, numel(peaksId2));
    for id = 1:2:numel(configuration.conditionEpochs)
        conditionBool(id:id + 1) = [configuration.conditionEpochs(id), ismember(allIds, time2id(t, configuration.conditionEpochs{id + 1}))];
        k = conditionBool{id + 1}(peaksId2);
        conditionIds(k) = (id + 1) / 2;
    end
    
    % Plot raw signal and bleaching.
    figure();
    ax.pre = axes('XTick', []);
    hold(ax.pre, 'all');
    ax.pre.Position = [ax.pre.Position(1), 0.7, ax.pre.Position(3), 0.3];
    plot(ax.pre, t, s);
    plot(ax.pre, t, sBleaching, '--');
    legend(ax.pre, {'Signal', 'Bleaching'});
    
    % Plot corrected signal, low-pass filtered signal and peaks.
    ax.post = axes('XTick', []);
    hold(ax.post, 'all');
    ax.post.Position = [ax.post.Position(1), 0.40, ax.post.Position(3), 0.25];
    plot(ax.post, t, z);
    plot(ax.post, t, signalLowPass, 'k');
    plot(ax.post, t([1, end]), threshold([1, 1]), 'k--');
    plot(ax.post, t(peaksId), z(peaksId), 'ro');
    legend(ax.post, {'z-score', 'low-pass', 'threshold', 'peaks'});
    
    % Plot traces around peaks.
    ax.dff = axes();
    xlabel(ax.dff, 'Time (s)');
    hold(ax.dff, 'all');
    ax.dff.Position = [ax.dff.Position(1), 0.1, ax.dff.Position(3), 0.25];
    t2 = [t(triggeredId) NaN(size(triggeredId, 1), 1)]';
    y2 = [dff(triggeredId) NaN(size(triggeredId, 1), 1)]';
    plot(ax.dff, t2(:), y2(:));
    plot(ax.dff, t(peaksId), dff(peaksId), 'ro');
    legend(ax.dff, {'df/f', 'peaks'});
    
    linkaxes([ax.pre, ax.post, ax.dff], 'x');
    axis(ax.pre, 'tight');
    
    % Plot triggered average.
    figure('name', 'FPA');
    ax.trigger = axes();
    hold(ax.trigger, 'all');
    uIds = unique(conditionIds);
    uIds = uIds(uIds > 0);
    triggerT = windowTemplateIds / configuration.resamplingFrequency;
    nPeaks = zeros(size(uIds));
    for i = 1:numel(uIds)
        id = uIds(i);
        nPeaks(i) = sum(conditionIds == id);
        data = dff(triggeredId(conditionIds == id, :));
        data = reshape(data, numel(data) / window, window);
        triggerMean = mean(data, 1);
        triggerStd = std(data, [], 1) / sqrt(size(data, 1));
        h1 = plot(triggerT, triggerMean);
        vertices = [triggerT; triggerMean + triggerStd];
        vertices = cat(2, vertices, [fliplr(triggerT); fliplr(triggerMean - triggerStd)])';
        faces = 1:2 * window;
        patch('Faces', faces, 'Vertices', vertices, 'FaceColor', h1.Color, 'EdgeColor', 'none', 'FaceAlpha', 0.10);
    end
    legend(gca, [configuration.conditionEpochs(2 * uIds - 1); arrayfun(@(s) sprintf('n=%i', s), nPeaks, 'UniformOutput', false)]);
    h = title('Triggered average for each condition \pm SEM');
    set(h, 'Interpreter', 'tex');
    xlabel('Time (s)');
    ylabel('df/f');
    axis(ax.trigger, 'tight');
end

function id = time2id(time_vector, time_limits)
    nEpochs = numel(time_limits) / 2;
    epochs = zeros(2, nEpochs);
    epochs(1:2:end) = arrayfun(@(l) find(time_vector >= l, 1, 'first'), time_limits(1:2:end));
    epochs(2:2:end) = arrayfun(@(l) find(time_vector <= l, 1, 'last'), time_limits(2:2:end));
    id = arrayfun(@(e) epochs(1, e):epochs(2, e), 1:nEpochs, 'UniformOutput', false);
    id = [id{:}];
end

function configuration = setDefault(configuration, fieldname, value)
    if ~isfield(configuration, fieldname)
        configuration.(fieldname) = value;
    end
end