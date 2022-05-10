% Add dependencies.
addpath('..');
addpath(genpath('../common'));

% Fiber-photometry data recorded with TDT DAQ.
inputFolder = 'C:\Users\Molina\Documents\public\data\HALO\FibrePhotometry\Gavin\FP_Jun2020/RA_FST_305-200318-133857';
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
configuration.conditionEpochs = {'Before', [626, 926], 'During1', [926, 1226]};
configuration.baselineEpochs = [100, 1500];
configuration.airPLS = false;
configuration.lowpassFrequency = 10;
configuration.peaksLowpassFrequency = 0.5;
configuration.baselineLowpassFrequency = 0.1;
configuration.resamplingFrequency = 100;
configuration.f0 = @median;
configuration.f1 = @mad;
configuration.threshold = {2.91, @mad, @median};
configuration.fitReference = true;
configuration.triggeredWindow = 10;

fpa = FPA(time, signal, reference, configuration);
cellfun(@warning, fpa.warnings);
fpa.plotTrace();