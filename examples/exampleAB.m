% Define epochs.
configuration = struct();
configuration.conditionEpochs = {'Baseline', [0, 800], 'Stroke', [800, 1600], 'Treatment', [1600, Inf]};

% DLC filename.
dlcFile = 'data/DLC.csv';
% Arena size (cm).
arenaSize = [100, 100];
% Number of pixels per cm.
ppcm = 10;
% Video frame rate.
framerate = 30;
% Minimum confidence level from DLC.
dlcThreshold = 0.90;
% Heat map saturation.
heatmapSaturation = 5;
heatmapBinCount = 20;

% Load behavioral data.
dlc = loadDLC(dlcFile);

% Choose x and y for body position.
x = dlc.bodypart1_x / ppcm;
y = dlc.bodypart1_y / ppcm;
p = dlc.bodypart1_p;
nFrames = size(dlc, 1);
time = colon(0, nFrames - 1)' / framerate;

% Filter.
cutoff = framerate / 2;
order = 9;
db = 20;
[b, a] = cheby2(order, db, cutoff / framerate / 2, 'low');
mask = p > dlcThreshold;
x2 = filtfilt(b, a, x(mask));
y2 = filtfilt(b, a, y(mask));
x2 = interp1(time(mask), x2, time);
y2 = interp1(time(mask), y2, time);

% Locomotion trace.
ax = NaN(1, 2);
ax(1) = subplot(1, 2, 1);
plot(x, y);
title('Raw');
xlabel('x (cm)');
ylabel('y (cm)');
axis('square');
ax(2) = subplot(1, 2, 2);
plot(x2, y2);
title('Processed');
xlabel('x (cm)');
axis('square');
linkaxes(ax);
axis(ax(1), 'tight');

% Locomotion over time.
delta = [0; sqrt(diff(x2) .^ 2 + diff(y2) .^ 2)];
figure()
plot(time, cumsum(delta));
axis('tight');
xlabel('time (s)');
ylabel('Locomotion (cm)');

% Total locomotion for individual epochs.
nEpochs = numel(configuration.conditionEpochs) / 2;
labels = configuration.conditionEpochs(1:2:end);
locomotion = zeros(nEpochs, 1);
for e = 1:nEpochs
    range = configuration.conditionEpochs{2 * e - 0};
    k = time2id(time, range);
    locomotion(e) = sum(delta(k));
end
figure()
bar(locomotion)
xticklabels(labels);
xlabel('Epoch');
ylabel('Locomotion (cm)');
xtickangle(45);
grid('on');
grid('minor');

% Locomotion heatmap.
xlims = [0, arenaSize(1)] + min(x);
ylims = [0, arenaSize(2)] + min(y);
xEdges = linspace(xlims(1), xlims(2), heatmapBinCount);
yEdges = linspace(ylims(1), ylims(2), heatmapBinCount);
[xGrid, yGrid] = meshgrid(xEdges, yEdges);
xGrid = xGrid(1:end - 1, 1:end - 1);
yGrid = yGrid(1:end - 1, 1:end - 1);
intensity = histcounts2(y2(:), x2(:), yEdges, xEdges);
figure();
pcolor(xGrid, yGrid, intensity);
shading('interp');
clim([0, prctile(intensity(:), 100 - heatmapSaturation)]);
axis('square');
xlabel('x (cm)');
ylabel('y (cm)');
