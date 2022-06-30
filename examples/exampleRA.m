% Fiber-photometry data recorded with TDT DAQ. (this extracts the data)
inputFolder = '../data/RA_Odor_604-201116-160100';
% Load 465nm and 405nm.
signalTitle = 'Dv1A';
referenceTitle = 'Dv2A';
data = loadTDT(inputFolder, {signalTitle, referenceTitle});
time = data(:, 1);
signal = data(:, 2);
reference = data(:, 3);

% Epoch definitions.
citralBaseline = [929 989];
citral = [989 1049 1049 1109 1109 1169 1169 1229 1229 1289 1825 1885];
citralPost = [1885 1945 1945 2005];
bobcatBaseline = [2145 2205];
bobcat =   [2205 2265 2265 2325 2325 2385 2385 2445 2445 2505 3028 3088];
bobcatPost = [3088 3148 3148 3208];

% Call FPA with given configuration.
configuration = struct();

configuration.f0 = {@mean, citralBaseline};
configuration.f1 = {@std, citralBaseline};
configuration.lowpassFrequency = 10;
configuration.resamplingFrequency = 100;
configuration.baselineEpochs = [0, 3500];
configuration.conditionEpochs = {'Citral baseline', citralBaseline, 'Citral', citral, 'Post citral', citralPost, 'Bobcat baseline', bobcatBaseline, 'Bobcat', bobcat, 'Post bobcat', bobcatPost};
configuration.peakWindow = [-5, 5];
configuration.events = [citral(1:2:end), bobcat(1:2:end)];
configuration.eventWindow = [-1, 29];
configuration.threshold = {2.91, @mad, @median};
configuration.fitReference = true;

% Call FPA with given configuration.
fpa = FPA(time, signal, reference, configuration);
cellfun(@warning, fpa.warnings);
fpa.plot();