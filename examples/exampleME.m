% Add dependencies.
addpath('..');
addpath(genpath('../common'));

% Fiber-photometry data recorded with Deisseroth's Multifiber.
inputDataFile = '../data/Multifiber.csv';
data = readtable(inputDataFile);
time = data.Time_410nm;
signal = data.MeanInt_470nm;
reference = data.MeanInt_410nm;
cleanEpochs = [20, 340];

configuration = struct();
configuration.conditionEpochs = {'Data', cleanEpochs};
configuration.baselineEpochs = cleanEpochs;
configuration.baselineLowpassFrequency = 0.1;
configuration.airPLS = true;
configuration.lowpassFrequency = 2;
configuration.peaksLowpassFrequency = 0.5;
configuration.threshold = {2.91, @mad, @median};
configuration.f0 = {@median, cleanEpochs};
configuration.f1 = {@std, cleanEpochs};

% Call FPA with given configuration.
fpa = FPA(time, signal, reference, configuration);
cellfun(@warning, fpa.warnings);

figure();
hold('on');
plot(fpa.time, fpa.signal, 'DisplayName', 'Corrected signal');
plot(fpa.time, fpa.reference, 'DisplayName', 'Corrected reference');
xlim(configuration.baselineEpochs);
xlabel('Time (s)');
legend('show');