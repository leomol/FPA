% Add dependencies.
addpath('..');
addpath(genpath('../common'));

% Fiber-photometry data recorded with TDT DAQ.
inputFolder = '../data/RA_FST_305-200318-133857';
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
configuration.threshold = {2.91, @mad, @median};
configuration.fitReference = true;
configuration.triggeredWindow = 10;
configuration.plot = true;

FPA(time, signal, reference, configuration);