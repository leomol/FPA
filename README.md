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
- Resample data to a target frequency.
- Correct for photo-bleaching with an exponential decay fit on the low-pass filtered signal.
- Normalize signal to reference.
- Compute df/f in a (moving) window.
- Band-pass filter data.
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