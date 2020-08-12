# Fiber-photometry analysis
MATLAB scripts to plot data from a fiber-photometry recording.

## Prerequisites
* [MATLAB][MATLAB] (last tested with R2020a)

## Installation
* Install MATLAB with the following toolboxes:
	* Curve Fitting Toolbox
	* Signal Processing Toolbox
	* Statistics and Machine Learning Toolbox
* Download and extract these scripts to Documents/MATLAB folder.

## Usage
For usage, see and edit examples according to your data recordings.

```matlab
FPA(time, signal, reference, configuration)
```

Correct signal from bleaching and artifacts, normalize and detect peaks of spontaneous activity based on the given parameters. Signal and reference are column vectors.

[![FPA demo](fpa-screenshot.png)](https://drive.google.com/file/d/1OXrwykbzTlqiQ13bCYg5v_xJNJIpqeb0)

## Analysis
Overall analysis steps:
- Resample to target frequency.
- Remove artifacts: Replace artifacts with lines in flagged regions.
- Correct for bleaching: Fit an exponential decay in low-pass data.
- Correct for motion artifacts: Subtract bleaching corrected signals.
- Compute df/f or z-score according to settings.
- Find peaks of spontaneous activity in band-pass signal.
- Plot 1:
  - Raw signal and bleaching fit.
  - Band-pass signal and peaks.
  - Motion and bleaching corrected signal.
- Plot 2:
  - Power spectrum of each epoch.
- Plot 3:
  - Boxplot.
- Plot 4:
  - Triggered average of spontaneous activity (if any peaks are found).

## Configuration

`configuration` is a struct with the following fields (defaults are used for missing fields):
- `conditionEpochs` - Epochs for different conditions: `{'epoch1', [start1, end1, start2, end2, ...], 'epoch2', ...}`
- `bleachingEpochs` - Time epochs (s) to include for bleaching correction.
- `artifactEpochs` - Time epochs (s) to remove.
- `resamplingFrequency` - Resampling frequency (Hz).
- `dffLowpassFrequency` - Lowpass frequency to filter `df/f`.
- `peaksBandpassFrequency` - Low/High frequencies to compute peaks.
- `bleachingLowpassFrequency` - Lowpass frequency to detect bleaching decay.
- `thresholdingFunction` - One of `@mad`, `@std`.
- `thresholdFactor` - Thresholding cut-off.
- `triggeredWindow` - Length of time to capture around each peak of spontaneous activity.
- `dffEpochs` - Time epochs (s) to include for `df/f` normalization.
- `f0Function` - One of `@movmean`, `@movmedian`, `@movmin`.
- `f0Window` - Length of the moving window to calculate `f0`.
- `f1Function` - One of `@movmean`, `@movmedian`, `@movmin`, `@movstd`.
- `f1Window` - Length of the moving window to calculate `f1`.

Peaks are calculated after bandpass filtering `df/f`. Everything else is calculated after lowpass filtering `df/f`.

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

### Example case 1:
Use the first minute as the baseline for the rest of the data.
```matlab
configuration.f0Function = @movmean;
configuration.f1Function = @movstd;
configuration.dffEpochs = [0, 60]
```
`f0Window` and `f1Window` are ignored because `dffEpochs` takes precedence.

### Example case 2:
Compute a moving baseline that is 1min wide.
```matlab
configuration.f0Function = @movmean;
configuration.f1Function = @movstd;
configuration.f0Window = 60;
configuration.f1Window = 60;
```
`dffEpochs` must be empty (`configuration.dffEpochs = [];`) which is the default.

### Example case 3:
Compute a common baseline from the entire recording.
```matlab
configuration.f0Function = @movmean;
configuration.f1Function = @movstd;
configuration.f0Window = Inf;
configuration.f1Window = Inf;
```
`dffEpochs` must be empty.

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

### Load an `abf` datafile acquired with Axon DAQ
```matlab
[data, units, names] = loadABF(filename)
```

### Load a data folder acquired with TDT DAQ
```matlab
[data, names] = loadTDT(folder)
```

Parses a project folder stored by a TDT DAQ and returns a `data` matrix where columns correspond to channels listed in `names`.

### Load a `mat` datafile acquired with LabChart
```matlab
[data, comments, units, titles] = loadLabChart(filename)
```

Loads a `mat` file exported by LabChart in an intuitive manner.
Numerical data and comments are organized in 2D structures where row indices select channels and column indices select blocks:
  - Numerical data:
```matlab
    data(channel, block).time
    data(channel, block).data
```
    
  - Comments:
```matlab
    comments(channel, block).time
    comments(channel, block).text
```

Units are organized in a cell array where row indices select channels and column indices select blocks:
   - Data units:
```matlab
    units{channel, block}
```
  
Channel titles are organized in a cell array:
```matlab
    titles{channel}
```

### Load TTL data logged by Inscopix DAQ
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
© 2018-2020 [Leonardo Molina][Leonardo Molina]

This project is licensed under the [GNU GPLv3 License][LICENSE.md].

[Leonardo Molina]: https://github.com/leomol
[MATLAB]: https://www.mathworks.com/downloads/
[LICENSE.md]: LICENSE.md