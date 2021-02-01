# Fiber-photometry analysis
MATLAB scripts to plot data from a fiber-photometry recording.

## Prerequisites
* [MATLAB][MATLAB] (last tested with R2020b)

## Installation
* Install MATLAB with the following toolboxes:
	* Curve Fitting Toolbox
	* Signal Processing Toolbox
	* Statistics and Machine Learning Toolbox
* Download and extract these scripts to Documents/MATLAB folder.

## Usage
See and edit examples according to your experimental setup.

### Script
```matlab
results = FPA(time, signal, reference, configuration);
```

### GUI
The GUI is under development and has limitted functionality compared to the script.

```matlab
GUI();
```

Remove baseline from data, correct from motion artifacts; normalize, filter, and detect peaks of spontaneous activity in user defined epochs.

## Analysis
Processing steps:
- Resample data to target frequency.
- Replace artifacts with linear interpolation in flagged regions.
- Baseline correction modeled as an exponential decay of the low-pass filtered data.
- Correct for motion artifacts by subtracting reference to signal, after a polynomial fit.
- Remove fast oscillations with a low-pass filter.
- Normalize data as df/f or z-score according to settings.
- Find peaks of spontaneous activity in low-pass filtered data.
- Plot 1:
  - Raw signal and reference, and baseline model.
  - Baseline corrected signal and reference.
  - Motion correction.
  - Normalization.
  - Peaks.
- Plot 2:
  - Power spectrum for each epoch.
- Plot 3:
  - Boxplot.
- Plot 4:
  - Triggered average of spontaneous activity (if any peaks are found).

## Configuration

`configuration` is a struct with the following fields (defaults are used for missing fields):
- `conditionEpochs` - Epochs for different conditions: `{'epoch1', [start1, end1, start2, end2, ...], 'epoch2', ...}`
- `baselineEpochs` - Time epochs (s) to include for baseline correction.
- `baselineLowpassFrequency` - Frequency representative of baseline.
- `airPLS` - Baseline correction for all data using airPLS (true, false, or airPLS inputs).
- `artifactEpochs` - Time epochs (s) to remove.
- `resamplingFrequency` - Resampling frequency (Hz).
- `lowpassFrequency` - Lowest frequency permitted in normalized signal.
- `peaksLowpassFrequency` - Lowest frequency to detect peaks.
- `thresholdingFunction` - `@mad`, `@std`, ...
- `thresholdFactor` - Thresholding cut-off.
- `triggeredWindow` - Length of time to capture around each peak of spontaneous activity.
- `fitReference` - Shift and scale reference to fit signal.

See source code for default values:
```matlab
edit('FPA');
```

Units for time and frequency are seconds and hertz respectively.

## Normalization recipes

Normalization is calculated as `(f - f0) / f1` where `f0` and `f1` can be data provided by the user or calculated using given functions:

### f0 and f1 are common to all datapoints and calculated from all data:
#### df/f:
```matlab
configuration.f0 = @mean;
configuration.f1 = @mean;
```

#### z-score:
```matlab
configuration.f0 = @mean;
configuration.f1 = @std;
```

#### z-score - alternative 1 (default):
```matlab
configuration.f0 = @median;
configuration.f1 = @mad;
```

#### z-score - alternative 2:
```matlab
configuration.f0 = @median;
configuration.f1 = @std;
```

### f0 and f1 are common to all data points and calculated at given epochs:
#### df/f:
```matlab
epochs = [0, 100, 500, 550, 1000, Inf]
configuration.f0 = {@mean, epochs};
configuration.f1 = {@mean, epochs};
```

### f0 and f1 are calculated for each data point based on a moving window:
#### df/f:
```matlab
window = 60;
configuration.f0 = {@movmean, window};
configuration.f1 = {@movmean, window};
```

[further combinations possible with `@movmean`, `@movmedian`, `@movstd`, `@movmad`, `@mov`...].

### Normalization from given data:
```matlab
f0 = ones(size(time));
f1 = ones(size(time)) * 10;
configuration.f0 = f0;
configuration.f1 = f1;
```

## Data loaders
### Load a CSV file (e.g. data acquired with Doric, Multifiber, or Inscopix DAQ)
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

### Load LabChart data
```matlab
[time, data, units, names, comments] = loadLabChart(filename)
```
Loads a `mat` file exported by LabChart in an intuitive manner.

```matlab
[time, data, units, names] = loadAdicht(filename)
```
Loads an `aditch` file.

Data are organized in 2D cells where row indices select channels and column indices select blocks:
```matlab
    time{channel, block}
    data{channel, block}
```

Units are organized in a cell array where row indices select channels and column indices select blocks:
   - Data units:
```matlab
    units{channel, block}
```
  
Channel names are organized in a cell array:
```matlab
    names{channel}
```

### Load TTL data logged by Inscopix DAQ
```matlab
ttl = loadInscopixTTL(filename)
```

Returns timestamps where pin `IO1` changes from low to high state.

### Load behavioral events detected by CleverSys
```matlab
[labels, start, finish, distance, speed] = loadCleverSys(filename, <sheet>)
epochs = loadCleverSys(filename, <sheet>)
```

Returns a list of event epochs generated by CleverSys.

## Version History
* 2021-01-15: Reference is fitted to signal. Band-pass replaced by low-pass. GUI example supports artifact removal.
* 2019-08-22: Library and example code. Added loaders for multiple data acquisition systems.
* 2019-02-20: Initial release.

## License
© 2018-2021 [Leonardo Molina][Leonardo Molina]

This project is licensed under the [GNU GPLv3 License][LICENSE.md].

[Leonardo Molina]: https://github.com/leomol
[MATLAB]: https://www.mathworks.com/downloads/
[LICENSE.md]: LICENSE.md