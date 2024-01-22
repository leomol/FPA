# Fiber photometry analysis
MATLAB toolbox to process, plot and export fiber photometry data.

## Prerequisites
- [MATLAB][MATLAB] (last tested with `R2023b`)

## Installation
- Install [MATLAB][MATLAB] with the following toolboxes:
	- Curve Fitting Toolbox (for `curvefit`, `fittype`, ...)
	- Signal Processing Toolbox (for `designfilt`, `filtfilt`, ...)
	- Statistics and Machine Learning Toolbox (for `zscore`, `mad`, ...)
- Download and extract the [FPA library][FPA] to the `Documents/MATLAB` folder

## Usage
- Run `startup.m` to add FPA dependencies to the MATLAB search path.
- Create a new script using the examples below as a template for your experimental setup, then run it.

### Example usage - minimal script
Call FPA using all default settings then call `plotTrace` (one of the functions for [plotting](#plotting-functions) or [exporting](#exporting-functions)).

Change the data loader (and data columns) according to your data acquisition system and data layout.

```MATLAB
filename = 'data.csv';
data = CSV.load(filename);
time = data(:, 1);
signal = data(:, 2);
reference = data(:, 3);

fpa = FPA(time, signal, reference);
fpa.plotTrace();
```

### Modifying default configuration
Configure each of the analysis steps. With the help of this documentation, you may adjust some values or whole functions altogether to customize the analysis.

```MATLAB
config = struct();
config.epochs = {'Data', [-Inf, Inf]};
config.resampleData = @(fpa) fpa.resample(100);
config.trimSignal = @(fpa, time, data) fpa.interpolate(time, data, []);
config.trimReference = @(fpa, time, data) fpa.interpolate(time, data, []);
config.smoothSignal = @(fpa, time, data) fpa.lowpass(time, data, 0.1);
config.smoothReference = @(fpa, time, data) fpa.lowpass(time, data, 0.1);
config.modelSignal = @(fpa, time, data) fpa.fit(time, data, 'exp1', [-Inf, Inf]);
config.modelReference = @(fpa, time, data) fpa.fit(time, data, 'exp1', [-Inf, Inf]);
config.correctSignal = @(fpa) fpa.signalTrimmed - fpa.signalModeled;
config.correctReference = @(fpa) fpa.referenceTrimmed - fpa.referenceModeled;
config.standardizeSignal = @(fpa) zscore(fpa.signalCorrected);
config.standardizeReference = @(fpa) zscore(fpa.referenceCorrected);
config.fitReference = @(fpa) fpa.fitReferenceToSignal(0.1, [-Inf, Inf]);
config.getF = @(fpa) fpa.signalStandardized - fpa.referenceFitted;
config.smoothF = @(fpa, time, data) fpa.lowpass(time, data, 10);
config.normalizeF = @(fpa, time, data) fpa.normalize(time, data, @median, @mad);
config.normalizeEvents = @(fpa, time, data) fpa.normalize(time, data, {@mean, [-Inf, 0]}, 1);
config.normalizePeaks = @(fpa, time, data) fpa.normalize(time, data, {@mean, [-Inf, 0]}, 1);
config.getPeaks = @(fpa) fpa.findPeaks('prominence', 0.5, median(fpa.fNormalized) + 2.91 * mad(fpa.fNormalized));

filename = 'data.csv';
data = CSV.load(filename);
time = data(:, 1);
signal = data(:, 2);
reference = data(:, 3);

fpa = FPA(time, signal, reference, config);
fpa.plotTrace();
```

## FPA inputs
```MATLAB
fpa = FPA(time, signal, reference);
```
or
```MATLAB
fpa = FPA(time, signal, reference, configuration);
``` 
By default, FPA removes the photobleaching effect, corrects from motion artifacts, normalizes, and filters data.

The analysis steps can be adjusted with the `configuration`.

# Analysis steps and results

The processing steps are functions applied to 405 and 465nm data to ultimately yield `fNormalized` (or "`df/f`").

The steps are defined in the `configuration` and results are saved as properties in the `FPA` object.

|        Step            |                       Result                             |
|:----------------------:|:--------------------------------------------------------:|
| `resampleData`         | `timeResampled`, `signalResampled`, `referenceResampled` |
| `trimSignal`           | `signalTrimmed`                                          |
| `trimReference`        | `referenceTrimmed`                                       |
| `smoothSignal`         | `signalSmoothed`                                         |
| `smoothReference`      | `referenceSmoothed`                                      |
| `modelSignal`          | `signalModeled`                                          |
| `modelReference`       | `referenceModeled`                                       |
| `correctSignal`        | `signalCorrected`                                        |
| `correctReference`     | `referenceCorrected`                                     |
| `standardizeSignal`    | `signalStandardized`                                     |
| `standardizeReference` | `referenceStandardized`                                  |
| `fitReference`         | `referenceFitted`                                        |
| `getF`                 | `f`                                                      |
| `smoothF`              | `fSmooth`                                                |
| `normalizeF`           | `fNormalized`                                            |
| `getPeaks`             | `peakIds`, `peakLabels` `peakCounts`                     |
|          *             | `duration`                                               |
|          *             | `area`                                                   |

\* Data available after all the processing steps have been executed.

The user continues the analysis by plotting and exporting data according to epoch definitions or event-triggers.

### Nomenclature
The following are frequently used as verbs or adjectives in function names, property names, file headers, plots, and the documentation:
  - `signal`: 465nm trace
  - `reference`: 405nm trace
  - `f`: Fluorescence trace deduced from operations on the signal and reference traces
  - `dff`: Normalized `f`, which may correspond to `df/f`, `z-score`, or whichever normalization is provided
  - `smooth`: Lowpass filter for photobleaching detection
  - `model`: Obtain the underlying photobleaching baseline
  - `trim`: Remove manually identified artifacts
  - `correct`: Subtract photobleaching model from data
  - `standardize` / `normalize`: Apply a normalization function to the data

Property names that modify data are written in imperative form whereas the resulting data has an adjective appended to it:
  - `smoothSignal` is a function that smooths the 465nm trace
  - `signalSmoothed` is the smoothed 465nm trace
  - `timeResampled` is time after resampling and is common to all processing steps after the `resampling` step

`time`, `signal`, `reference` are raw data, prior to resampling, hence with a number of data points incompatible with any other processing step.

While data undergoes multiple processing steps, the adjectives attached to them only reflect the last operation performed. If such step were disabled in the pipeline, the analysis would still contain data from such step but the values would be identical to the previous step. Normally,

  `signal >> signalResampled >> signalTrimmed >> signalCorrected >> signalStandardized`,

  If `trimSignal` was disabled, `signalTrimmed` would contain the same data as `signalResampled`.

## Configuration
A default behavior applies for any function or parameter not defined in the `configuration`.

FPA helper functions used in the processing steps can be listed with `methods(FPA)` or by executing `FPA.defaults`

To disable any processing step, assign `[]` to the processing function.

Note that input and output units for time and frequency are given in seconds and hertz respectively.

The default processing steps are listed below.

### `resampleData`
Function that returns resampled `time`, `signal`, and `reference`.

The result is saved into `fpa.timeResampled`, `fpa.signalResampled`, and `fpa.referenceResampled`.

#### Default
Resample data to `100Hz` using a polyphase antialiasing filter. The request is ignored if the sampling rate is lower than this frequency.

#### Customization example
Resample data to 20 Hz reusing FPA's resampling method:
```MATLAB
config.resampleData = @(fpa) fpa.resample(20);
```

### `trimSignal`
Function that removes artifacts from the signal.

The result is saved into `fpa.signalTrimmed`.

#### Default
Ignored.

#### Customization example
Remove visually detected noise by applying a linear interpolation from time 0 to 1 then from 100 to 110 seconds over `fpa.signalResampled`:
```MATLAB
config.trimSignal = @(fpa, time, data) fpa.interpolate(time, data, [0, 1, 100, 110]);
```

### `trimReference`
Same as above.

### `smoothSignal`
Function that smooths `fpa.signalTrimmed`, for the purposes of modeling photobleaching.

The result is saved it into `fpa.signalSmoothed`.

#### Default
Apply a 12th order, lowpass, butter, digital filter with half power frequency of `0.1Hz` to the signal.

#### Example customization
Lowpass filter the signal to `0.5Hz`:
```MATLAB
config.smoothSignal = @(fpa, time, data) fpa.lowpass(time, data, 0.5);
```

### `smoothReference`
Same as above.

### `modelSignal`
Function that models photobleaching by fitting a curve to `fpa.signalSmoothed`.

The result is saved into `fpa.signalModeled`.

#### Default
Fit an exponential decay to the whole signal.

#### Example customization
Apply an exponential decay using data from the minute 0 to 10 and 50 to 60:
```MATLAB
config.modelSignal = @(fpa, time, data) fpa.fit(time, data, 'exp1', [0, 10, 50 60] * 60);
```

Using a linear fit instead:
```MATLAB
config.modelSignal = @(fpa, time, data) fpa.fit(time, data, 'poly1', [0, 10, 50 60] * 60);
```

### `modelReference`
Same as above.

### `correctSignal`
Function that removes photobleaching from the signal.

The result is saved into `fpa.signalCorrected`.

#### Default
Subtract `fpa.signalModeled` from `fpa.signalTrimmed`.

#### Example customization
Subtract the mean value calculated from a different recording.
```MATLAB
config.correctSignal = @(fpa) fpa.signalModeled - mean(data);
```
### `correctReference`
Same as above.

### `standardizeSignal`
Function that standardizes the signal.

The result is saved into `fpa.signalStandardized`

#### Default
Apply a `z-score` function to `fpa.signalCorrected`

#### Example customization
Apply a modified `z-score` function instead.
```MATLAB
config.standardizeSignal = @(fpa) (fpa.signalCorrected - median(fpa.signalCorrected)) ./ mad(fpa.signalCorrected);
```

### `standardizeReference`
Same as above

### `fitReference`
Function that rescales the reference so that it is comparable to the signal.

The result is saved into `fpa.referenceFitted`. 

#### Default
Fits a 0.1 Hz lowpass version of `fpa.referenceStandardize` to `fpa.signalStandardize` using a bisquare, linear regression over the whole range.

#### Example customization
Apply a non-negative, least-squares fit at the given epochs, after lowpass filtering to 0.5 Hz:
```MATLAB
config.fitReference = @(fpa) fpa.fitReferenceToSignal(0.5, [0, 300, 1200, 1500]);
```

### `getF`
Function that returns the motion- and bleaching corrected data.

The result is saved into `fpa.f`

#### Default
Subtract `fpa.referenceFitted` from `fpa.signalCorrected`.

#### Example customization
```MATLAB
config.getF = @(fpa) fpa.signalCorrected - fpa.referenceFitted;
```

### `smoothF`
Function that smooths `f`.

The result is saved into `fpa.fSmoothed`.

#### Default
Smooth `fpa.f` using a 12th order, lowpass, butter, digital filter with half power frequency of `10Hz`.

#### Example customization.
Smooth `f` using a lowpass filter with a frequency of `15Hz`:
```MATLAB
config.smoothF = @(fpa, time, data) fpa.lowpass(time, data, 15);
```

### `normalizeF`
Function that normalizes `f`. See the [Normalization section](#Normalization).

#### Default
Apply a modified `z-score`

The result is saved into `fpa.fNormalized`

#### Example customization
Apply a modified `z-score`:
```MATLAB
config.normalizeF = @(fpa, time, data) fpa.normalize(time, data, @median, @mad);
```

### `normalizeEvents`
Function that normalizes each event trace. See `config.normalizeF` above and the [Normalization section](#Normalization).

### `normalizePeaks`
Function that normalizes each peak trace. See `config.normalizeF` above and the [Normalization section](#Normalization).

#### Default
Subtract the mean value from the event trigger time (from `-Inf` to `0`).
```MATLAB
config.normalizePeaks = @(fpa, time, data) fpa.normalize(time, data, {@mean, [-Inf, 0]}, 1);
```

### `getPeaks`
Function to detect peaks for peak-triggered averages.

#### Default
Find peaks using the `prominence` option with at least 0.5 second separation, and a threshold of 2.91 median absolute deviations from the median
```MATLAB
config.getPeaks = @(fpa) fpa.findPeaks('prominence', 0.5, median(fpa.fNormalized) + 2.91 * mad(fpa.fNormalized));
```

#### Example customization
Find peaks using the `height` option with at least 1.0 second separation, and a threshold of 2 standard deviations from the mean
```MATLAB
config.getPeaks = @(fpa) fpa.findPeaks('height', 1.0, mean(fpa.fNormalized) + 2 * std(fpa.fNormalized));
```

### `epochs`
Structure with epoch definitions according to the experimental design.
Plots and output data will be parsed out and averaged according to these epochs.

#### Default
A single epoch called `Data` spanning the whole recording.
```MATLAB
config.epochs = struct('Data', [-Inf, Inf]);
```

#### Example customization
- Define a 10-minute habituation period at the start (0 to 600).
- Define 2x 5 second air puffs separated by 30 seconds at minute 10 (600 to 605 and 630 to 635).
- Define an injection period at minute 15, lasting 10 seconds (900 to 910).
- Define 2x 5 second air puffs separated by 30 seconds at minute 20 (1200 to 1205 and 1230 to 1235).

```MATLAB
config.epochs = struct('Habituation', [0, 600], 'AirPuffPre', [600, 605, 630, 635], 'Injection', [900, 910], 'AirPuffPost', [1200, 1205, 1230, 1235]);
```

# Plotting functions
- `plotAUC([normalized])` - Plot area under the curve which corresponds to the mean value when `normalized` is set to `true` (default).

- `plotPowerSpectrum([overlap], [window])` - Plot power spectrum via Welch's method, using a `window` of a given size (defaults to 10s). Set `overlap` to `true` to plot different epochs on top of each other (default) or `false` to separate in subplots.

- `plotStatistics()` - Make a boxplot per epoch.

- `plotTrace()` - Plot most analysis steps:
  - Raw data and photobleaching model
  - Photobleaching corrected data
  - Standardized data
  - Motion corrected data
  - Normalized data and peaks detected

- `plotEvents(eventTimes, [window])` - Plot event-triggered average for time points given by `eventTimes`, using a `window` of a given size (defaults to 5s).

- `plotEventHeatmap(eventTimes, [window])` - Similar to `plotEvents` but displaying each trigger using a heatmap plot.

- `plotPeaks([window])` - Plot peak-triggered average for peak times detected according to `config.getPeaks`, using a `window` of a given size (defaults to 5s).

- `plotPeakHeatmap([window])` - Similar to `plotPeaks` but displaying each trigger using a heatmap plot.

# Exporting functions
- `exportF(filename)` - Export time vs f  (both normalized and unnormalized).

- `exportStatistics(filename)` - Export AUC, sum, mean, median, min and max for normalized f and epoch duration and peak counts.

- `exportEvents(filename, eventTimes, [window], [asColumns])` - Export event-triggered traces with corresponding epoch name for time points given by `eventTimes`, using a `window` of a given size (defaults to 5s). When `asColumns` is set to `true` (default), the output `csv` file will show data in columns, otherwise as rows.

- `exportEventAverage(filename, eventTimes, [window], [asColumns])` - Similar to `exportEvents` but traces of the same epoch are averaged together.

- `exportPeaks(filename, [window], [asColumns])` - Export spike-triggered traces with corresponding epoch name for peak times detected according to `config.getPeaks`, using a `window` of a given size (defaults to 5s). When `asColumns` is set to `true` (default), the output `csv` file will show data in columns, otherwise as rows.

- `exportPeakAverage(filename, [window], [asColumns])` - Similar to `exportPeaks` but traces of the same epoch are averaged together.

# Helper functions
FPA methods to customize the configuration:
- `findPeaks(type, separation, peakThreshold)` - Find peaks in data.
- `fit(time, data, fitType, epochs)` - Fit curve to trace of choice according to the fitting type.
- `fitReferenceToSignal(targetFrequency, epochs)` - FFit reference to signal using non-negative, least-squares fit, at the given epochs.
- `get(traceName, epochs)` - Get trace according to data choice and epoch range.
- `ids(epochs)` - Get indices relative to `fpa.timeResampled` for the given epochs.
- `interpolate(time, data, epochs)` - Interpolate trace linearly according to data choice and epoch range.
- `lowpass(time, data, frequency)` - Low-pass filter trace according to data choice and frequency.
- `normalize(time, data, f0, f1)` - Normalize data according to parameters f0 and f1.
- `resample(frequency)` - Resample signal and reference according to frequency.

# Normalization
Normalization is calculated as `(f - f0) / f1` where `f0` and `f1` are parameters provided by the user.

`f0` and `f1` can be one of the following:
- A numeric value. For example, `0` or `1`.
- An array of numeric values of the same size of `f`. For example, `[0, 0, 0, ...]`
- A `function` to apply to the whole trace to obtain a single numeric value. For example, `@mean`, `@median`, `@std`, `@mad`
- A `function` together with time epochs to obtain a single numeric value from such epochs. For example, `{@mean, [start1, stop1, start2, end2, ...]}`
- A _moving_ function such as `@movmean`, `@movstd`, ... together with a window size to obtain a numeric value for each value of `f`.

## Normalization recipes
The following applies to `config.normalizeF`, `config.normalizePeaks`, and `config.normalizeEvents`

Apply a modified `z-score` where `f0` and `f1` are common to all datapoints and calculated from all data:
```MATLAB
@(fpa, time, data) fpa.normalize(time, data, @median, @mad);
```

Apply a regular `z-score`, calculate the mean from the first 10 minutes only:
```MATLAB
@(fpa, time, data) fpa.normalize(time, data, {@mean, [0, 600]}, @std);
```

Apply `df/f`:
```MATLAB
@(fpa, time, data) fpa.normalize(time, data, @mean, @mean);
```

Apply `df/f` using a moving window of 1 minute:
```MATLAB
@(fpa, time, data) fpa.normalize(time, data, {@movmean, 60}, {@movmean, 60});
```

further combinations possible with `@movmean`, `@movmedian`, `@movstd`, `@movmad`, ...

Normalize using known values
```MATLAB
@(fpa, time, data) (data – value1) / value2;
@(fpa, time, data) (data – array1) / array2;
```



## Data loaders
### Load a CSV file (e.g. data acquired with Doric, Multifiber, or Inscopix DAQ)
```MATLAB
[data, names] = CSV.load(filename)
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
```MATLAB
[data, names, sheetName] = XLS.load(filename, <sheetNumber|sheetName>)
```

Returns a `data` matrix where columns correspond to channels listed in `names`. A sheet number or a sheet name can be passed as the second parameter (default is the first sheet in the document). The sheet name is also returned.

Expected format for each sheed in a document:

| time | name 1 | name 2 | ... | name N |
|:----:|:-----: |:------:|:---:|:------:|
|  t1  |   a1   |   b2   | ... |   zN   |
|  ... |   ...  |  ....  | ... |   ...  |
|  tM  |   aM   |   bM   | ... |   zM   |

### Load Doric's HDF-formatted files (`.doric`)
```MATLAB
data = Doric.load(filename, [datasets])
```

Returns a `data` matrix where columns correspond to the `datasets` provided. When `datasets` is not provided, the function attempts to retrieve data assuming the typical layout.

Available `datasets` can be retrieved using `DoricStudio` or using FPA's function `Doric.getDatasets(filename)`, from such list, you can identify the name of the required datasets. Example `dataset` names:

```
datasets = {
  '/DataAcquisition/FPConsole/Signals/Series0001/AIN01xAOUT01-LockIn/Time'
  '/DataAcquisition/FPConsole/Signals/Series0001/AIN01xAOUT01-LockIn/Values'
  '/DataAcquisition/FPConsole/Signals/Series0001/AIN02xAOUT02-LockIn/Values'
};
```
or
```
datasets = {
  '/Traces/Console/Time(s)/Console_time(s)'
  '/Traces/Console/AIn-1 - Dem (AOut-1)/AIn-1 - Dem (AOut-1)'
  '/Traces/Console/AIn-2 - Dem (AOut-2)/AIn-2 - Dem (AOut-2)'
};
```

### Load an `abf` datafile acquired with Axon DAQ
```MATLAB
[data, units, names] = ABF.load(filename)
```

### Load a data folder acquired with TDT DAQ
```MATLAB
[data, names] = TDT.load(folder)
```

Parses a project folder stored by a TDT DAQ and returns a `data` matrix where columns correspond to channels listed in `names`.

### Load LabChart data
```MATLAB
[time, data, units, names, comments] = LabChart.Mat.load(filename)
```
Loads a `mat` file exported by LabChart in an intuitive manner.

```MATLAB
[time, data, units, names] = LabChart.Adicht.load(filename)
```
Loads an `aditch` file.

Data are organized in 2D cells where row indices select channels and column indices select blocks:
```MATLAB
    time{channel, block}
    data{channel, block}
```

Units are organized in a cell array where row indices select channels and column indices select blocks:
   - Data units:
```MATLAB
    units{channel, block}
```
  
Channel names are organized in a cell array:
```MATLAB
    names{channel}
```

### Load TTL data logged by Inscopix DAQ
```MATLAB
ttl = Inscopix.loadTTL(filename)
```

Returns timestamps where pin `IO1` changes from low to high state.

### Load behavioral events detected by CleverSys
```MATLAB
[labels, start, stop, distance, speed] = CleverSys.load(filename, <sheet>)
epochs = CleverSys.load(filename, <sheet>)
```

Returns a list of event epochs generated by CleverSys.

### Load behavioral events detected by BORIS (aggregated and tabulated export types)
```MATLAB
[labels, start, stop] = Boris.load(filename)
epochs = Boris.load(filename)
```

Returns a list of event epochs scored with BORIS.


### Load DeepLabCut data
```MATLAB
data = DLC.load(filename)
```

Reads a DeepLabCut-generated position file in `csv` format into a MATLAB table.


## Changelog
See [Changelog](CHANGELOG.md)

## License
© 2018 [Leonardo Molina][Leonardo Molina]

This project is licensed under the [GNU GPLv3 License][LICENSE.md].

[Leonardo Molina]: https://github.com/leomol
[FPA]: https://github.com/leomol/FPA
[MATLAB]: https://www.mathworks.com/downloads/
[LICENSE.md]: LICENSE.md