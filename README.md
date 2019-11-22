# Fiber-photometry analysis
MATLAB scripts to plot data from a fiber-photometry recording.

## Prerequisites
* [MATLAB][MATLAB] (last tested with R2019a)

## Installation
* Install MATLAB with the following toolboxes:
	* Curve Fitting Toolbox
	* Signal Processing Toolbox
	* Statistics and Machine Learning Toolbox
* Download and extract these scripts to Documents/MATLAB folder.

## Usage
See examples/description/video below and/or edit `FPAexamples.m` according to your data recordings.

```matlab
FPA(time, signal, reference, configuration)
```

Correct signal from bleaching and artifacts, normalize and detect peaks of spontaneous activity based on the given parameters. Signal and reference are column vectors.

[![FPA demo](fpa-screenshot.png)](https://drive.google.com/file/d/1OXrwykbzTlqiQ13bCYg5v_xJNJIpqeb0)

## Analysis
Overall analysis steps:
- Resample data.
- Fit an exponential decay in a portion of the data representative of bleaching.
- Normalize signal to reference.
- Compute df/f in a (moving) window.
- Find peaks of spontaneous activity.
- Plot corrected signal and peaks; highlight epochs.
- Plot triggered averages.
- Plot power spectrum.
- Plot stats.

## Configuration

`configuration` is a struct with the following fields (defaults are used for missing fields):
- `conditionEpochs` - Epochs for different conditions: `{'epoch1', [start1, end1, start2, end2, ...], 'epoch2', ...}`
- `baselineEpochs` - Time epochs (s) to include for df/f normalization.
- `bleachingEpochs` - Time epochs (s) to include for bleaching correction.
- `artifactEpochs` - Time epochs (s) to remove.
- `resamplingFrequency` - Resampling frequency (Hz).
- `dffLowpassFrequency` - Lowpass frequency to filter df/f.
- `peaksBandpassFrequency` - Low/High frequencies to compute peaks.
- `bleachingLowpassFrequency` - Lowpass frequency to detect bleaching decay.
- `f0Function` - One of `@movmean`, `@movmedian`, `@movmin`.
- `f0Window` - Length of the moving window to calculate f0.
- `f1Function` - One of `@movmean`, `@movmedian`, `@movmin`, `@movstd`.
- `f1Window` - Length of the moving window to calculate f1.
- `thresholdingFunction` - One of `@mad`, `@std`.
- `thresholdFactor` - Thresholding cut-off.
- `triggeredWindow` - Length of time to capture around each peak of spontaneous activity.

Peaks are calculated after bandpass filtering df/f. Everything else is calculated after lowpass filtering df/f.

See source code for default values:
```matlab
edit('FPA');
```

Units for time and frequency are seconds and hertz respectively.

## Normalization recipes

`df/f` is calculated as `(f - f0) / f1`, where `f0` and `f1` change according to the configuration.

For example, you may want to see changes relative to a 10s moving window:
```matlab
configuration.f0Function = @movmean;
configuration.f1Function = @movmean;
configuration.f0Window = 10;
configuration.f1Window = 10;
```

... or you may want to see changes from the standard deviation of the whole recording:
```matlab
configuration.f0Function = @movmean;
configuration.f1Function = @movstd;
configuration.f0Window = Inf;
configuration.f1Window = Inf;
```

If baselineEpochs covers the entire data set (e.g. `[-Inf, Inf]`), `df/f` is calculated using a moving window. If baselineEpochs covers a portion of the data set (e.g. `[10, 100]`), `df/f` is calculated using a single window in such period.

## Examples

### Example 1 - Fiber-photometry data recorded with Doric DAQ
```matlab
inputDataFile = 'data/Doric.csv';
% Names of columns corresponding to 465nm and 405nm.
signalTitle = 'AIn-1 - Demodulated(Lock-In)';
referenceTitle = 'AIn-2 - Demodulated(Lock-In)';
% Configuration (help FPA).
configuration = struct();
configuration.conditionEpochs = {'Pre', [100, 220], 'During', [650, 890], 'Post', [1480, 1600]};
configuration.bleachingEpochs = [-Inf, 600, 960, Inf];
configuration.resamplingFrequency = 100;
configuration.f0Function = @movmean;
configuration.f0Window = 600;
configuration.f1Function = @movmean;
configuration.f1Window = 600;
configuration.thresholdingFunction = @mad;
configuration.thresholdFactor = 0.1;
% Load data (help loadData).
[data, names] = loadData(inputDataFile);
% Identify columns in data.
s = ismember(names, signalTitle);
r = ismember(names, referenceTitle);
time = data(:, 1);
signal = data(:, s);
reference = data(:, r);
% Call FPA with given configuration.
FPA(time, signal, reference, configuration);
```

### Example 2 - Fiber-photometry data with stimuli recorded with Inscopix
```matlab
inputDataFile = 'data/Inscopix.csv';
inputEventFile = 'data/InscopixTTL.csv';
configuration = struct();
configuration.f0Function = @movmean;
configuration.f0Window = 600;
configuration.f1Function = @movmean;
configuration.f1Window = 600;
configuration.thresholdingFunction = @mad;
configuration.thresholdFactor = 0.1;
configuration.triggeredWindow = 1;
configuration.dffLowpassFrequency = 2;
configuration.peaksBandpassFrequency = [0.2, 2];
data = loadData(inputDataFile);
time = data(:, 1);
signal = data(:, 2);
% Extract epochs from TTL output.
eventDuration = 3;
baselineOffset = -4;
events = loadInscopixTTL(inputEventFile);
events = [events; events + eventDuration];
epochs = {'Baseline', events + baselineOffset, 'Stimulation', events};
configuration.conditionEpochs = epochs;
FPA(time, signal, [], configuration);
```

### Example 3 - Fiber-photometry data recorded with Doric DAQ and behavioral data recorded with CleverSys

```matlab
% Fibre photometry recording file.
inputDataFile = 'data/Doric.csv';
% CleverSys event file in seconds and the name of the target sheet within.
inputEventFile = {'data/CleverSys.xlsx', 'Trial 1'};
signalTitle = 'AIn-1 - Demodulated(Lock-In)';
referenceTitle = 'AIn-2 - Demodulated(Lock-In)';
configuration = struct();
configuration.resamplingFrequency = 20;
configuration.f0Function = @movmean;
configuration.f0Window = 60;
configuration.f1Function = @movmean;
configuration.f1Window = 60;
configuration.baselineEpochs = [-Inf, 600];
configuration.artifactEpochs = [603, 620, 910, 915];
configuration.bleachingEpochs = [-Inf, 600, 960, Inf];
configuration.thresholdingFunction = @mad;
configuration.thresholdFactor = 0.1;
% Extract epochs from CleverSys output.
events = loadCleverSysEvents(inputEventFile{:});
eventNames = events.keys;
configuration.conditionEpochs = cellfun(@(eventName) {eventName, reshape([events(eventName).start, events(eventName).start + events(eventName).duration]', 1, 2 * numel(events(eventName).start))}, eventNames, 'UniformOutput', false);
configuration.conditionEpochs = cat(2, configuration.conditionEpochs{:});
[data, names] = loadData(inputDataFile);
s = ismember(names, signalTitle);
r = ismember(names, referenceTitle);
time = data(:, 1);
signal = data(:, s);
reference = data(:, r);
results = FPA(time, signal, reference, configuration);
% Save peak times to file.
[folder, basename] = fileparts(inputDataFile);
output = fullfile(folder, sprintf('%s peak-time.csv', basename));
fid = fopen(output, 'w');
fprintf(fid, 'Peak Time (s)\n');
fprintf(fid, '%.3f\n', results.time(results.peaksId));
fclose(fid);
```

### Example 4 - Fiber-photometry data recorded with TDT DAQ
```matlab
inputFolder = 'data/GP_PVN_13a-190531-122516';
signalTitle = 'Dv1A';
referenceTitle = 'Dv2A';
configuration = struct();
configuration.conditionEpochs = {'Baseline', [1, 900], 'Test', [1102, 1702]};
configuration.baselineEpochs = [-Inf, Inf];
configuration.bleachingEpochs = [1, 1748];
configuration.resamplingFrequency = 20;
configuration.f0Function = @movmean;
configuration.f0Window = 10;
configuration.f1Function = @movmean;
configuration.f1Window = 10;
configuration.thresholdingFunction = @mad;
configuration.thresholdFactor = 0.10;
data = loadTDT(inputFolder, {signalTitle, referenceTitle});
time = data(:, 1);
signal = data(:, 2);
reference = data(:, 3);
FPA(time, signal, reference, configuration);
```

## Data loaders
### Load a CSV file (e.g. data acquired with Doric or Inscopix DAQ)
```matlab
[data, names] = loadData(filename)
```

Returns a `data` matrix where columns correspond to channels listed in `names`, that is, time followed by other data columns.

Expected format #1:

| time | name 1 | name 2 | ... | name N |
|:----:|:-----: |:------:|:---:|:------:|
|  t1  |   a1   |   b2   | ... |   zN   |
|  ... |   ...  |  ....  | ... |   ...  |
|  tM  |   aM   |   bM   | ... |   zM   |

Expected format #2:

|      | name 1 | name 2 | ... | name N |
|:----:|:-----: |:------:|:---:|:------:|
| time | accepted\|rejected | accepted\|rejected | ... | accepted\|rejected |
|  t1  |   a1   |   b2   | ... |   zN   |
|  ... |   ...  |  ....  | ... |   ...  |
|  tM  |   aM   |   bM   | ... |   zM   |

### Load an XLS or XLSX file (e.g. data acquired with Doric or Inscopix DAQ)
```matlab
[data, names, sheetName] = loadData(filename, <sheetNumber|sheetName>)
```

Returns a `data` matrix where columns correspond to channels listed in `names`. A sheet number or a sheet name can be passed as the second parameter (default is the first sheet in the document). The sheet name is also returned.

Expected format for each sheed in a document:

| time | name 1 | name 2 | ... | name N |
|:----:|:-----: |:------:|:---:|:------:|
|  t1  |   a1   |   b2   | ... |   zN   |
|  ... |   ...  |  ....  | ... |   ...  |
|  tM  |   aM   |   bM   | ... |   zM   |

### Load a data folder acquired with TDT DAQ
```matlab
[data, names] = loadTDT(folder)
```

Parses a project folder stored by a TDT DAQ and returns a `data` matrix where columns correspond to channels listed in `names`.

### Load TTL logged by Inscopix DAQ
```matlab
ttl = loadInscopixTTL(filename)
```

Returns timestamps where pin `IO1` changes from low to high state.

### Load behavioral events detected by CleverSys
```matlab
events = loadCleverSysEvents(filename, sheet)
```

Returns a map of events generated by CleverSys.

## Version History
* 0.1.1: Library and example code. Added loaders for multiple data acquisition systems.
* 0.1.0: Initial release.

## License
© 2019 [Leonardo Molina][Leonardo Molina]

This project is licensed under the [GNU GPLv3 License][LICENSE.md].

[Leonardo Molina]: https://github.com/leomol
[MATLAB]: https://www.mathworks.com/downloads/
[LICENSE.md]: LICENSE.md