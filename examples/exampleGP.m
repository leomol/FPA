% Add dependencies.
addpath('..');
addpath(genpath('../common'));

% Fiber-photometry data recorded with TDT DAQ.
inputFolder = '../data/RA_FST_305-200318-133857/';
signalTitle = 'Dv1A';
referenceTitle = 'Dv2A';
data = loadTDT(inputFolder, {signalTitle, referenceTitle});
time = data(:, 1);
signal = data(:, 2);
reference = data(:, 3);

% Call FPA with given configuration.
configuration = struct();
configuration.f0 = @median;
configuration.f1 = @mad;
configuration.conditionEpochs = {'Before', [626, 926] 'During1', [926, 1226] 'During2', [1226, 1499] 'During3', [1499, 1799] 'After', [1799, 2099],'SBefore', [1806, 2106] 'SDuring1', [2106, 2406] 'SDuring2', [2406, 2700] 'SDuring3', [2700, 3000] 'SAfter', [3000, 3300]};
configuration.baselineEpochs = [100, 1500];
configuration.airPLS = false;
configuration.lowpassFrequency = 10;
configuration.peakslowpassFrequency = 0.5;
configuration.baselineLowpassFrequency = 0.1;
configuration.resamplingFrequency = 100;
configuration.f0 = @median;
configuration.f1 = @mad;
configuration.thresholdingFunction = @mad;
configuration.thresholdFactor = 2.91;
configuration.fitReference = true;
configuration.triggeredWindow = 10;
configuration.plot = true;

results = FPA(time, signal, reference, configuration);

% Area under the curve.
figure();
nEpochs = numel(configuration.conditionEpochs) / 2;
epochNames = configuration.conditionEpochs(1:2:2 * nEpochs);
area = zeros(1, nEpochs);
duration = zeros(1, nEpochs);
for e = 1:nEpochs
    ids = time2id(results.time, configuration.conditionEpochs{2 * e});
    dff = results.dff(ids);
    area(e) = sum(dff);
    duration(e) = numel(ids);
end
bar(1:nEpochs, area);
h = gca();
h.XTickLabel = epochNames;
title('dff/f - AUC');

figure();
bar(1:nEpochs, area ./ duration);
h = gca();
h.XTickLabel = configuration.conditionEpochs(1:2:2 * nEpochs);
title('dff/f - normalized AUC');