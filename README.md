# Fiber photometry analysis
MATLAB scripts to plot data from a fiber photometry recording.

## Prerequisites
* [MATLAB][MATLAB] (last tested with R2018a)

## Installation
* Install MATLAB.
* Download and extract these scripts to Documents/MATLAB folder.

## Usage
```matlab
	FPA(inputFile, configuration)
```
Plot spontaneous signal from a fiber-photometry experiment based on the values of the `configuration`.

`inputFile` is a CSV file containing time, signal and reference columns.

<!---
[![FPA demo](fpa-snapshop.png)](https://drive.google.com/file/d/19h34s5LPmWgZJFF17zxef8f8A4bYAu90)
-->

## Analysis
- Load and resample `inputFile`.
- Low-pass filter an artifact-free portion of the data and fit an exponential decay to correct for bleaching.
- Correct for movement artifacts with reference signal.
- Compute z-score over the whole recording and low-pass filter to detect peaks.
- Compute df/f in a moving time window to normalize traces around peaks.
- Compute triggered averages of spontaneous activity grouped by condition/epochs definition.

`configuration` is a struct with the following fields (defaults are used for missing fields):
- `timeTitle` - Exact title name of time column.
- `signalTitle` - Exact title name of time column.
- `referenceTitle` - Exact title name of time column.
- `resamplingFrequency` - Resampling frequency (Hz).
- `bleachingCorrectionEpochs` - Time epochs (s) to include for bleaching correction.
- `f0Function` - One of `@movmean`, `@movmedian`, `@movmin`.
- `f0Window` - Length of the moving window to calculate f0.
- `f1Function` - One of `@movmean`, `@movmedian`, `@movmin`, `@movstd`.
- `f1Window` - Length of the moving window to calculate f1.
- `peaksLowpassFrequency` - Frequency of lowpass filter to detect peaks.
- `thresholdingFunction` - One of `@mad`, `@std`.
- `thresholdFactor` - Thresholding cut-off.
- `conditionEpochs` - Epochs that involving different conditions: `{'epoch1', [start1, end1, start2, end2, ...], 'epoch2', ...}`
- `triggeredWindow` - Length of the time window around each peak.

See source code for default values.

Units for time and frequency are seconds and hertz respectively.

Normalization:
    `df/f` is calculated as `(f - f0) / f1`, where `f0` and `f1` are computed for each time point using function over a moving window with the given size (s).
	
    Recipe for `df/f`:
	```matlab
        configuration.f0Function = @movmean;
        configuration.f1Function = @movmean;
	```
    Recipe for z-score:
	```matlab
        configuration.f0Function = @movmean;
        configuration.f1Function = @movstd;
	```

## Example 1:
```matlab
    inputFile = 'data/tf_CP_m307_d2_2.csv';
    FPA(inputFile);
```

## Example 2:
```matlab
	inputFile = 'data/tf_CP_m307_d2_2.csv';
	configuration.timeTitle = 'Time(s)';
	configuration.signalTitle = 'AIn-2 - Demodulated(Lock-In)';
	configuration.referenceTitle = 'AIn-1 - Demodulated(Lock-In)';
	configuration.resamplingFrequency = 50;
	configuration.bleachingCorrectionEpochs = [-Inf, 600, 960, Inf];
	configuration.conditionEpochs = {'Condition A', [-Inf, 600], 'Condition B', [650, Inf]};
	configuration.triggeredWindow = 10;
	configuration.f0Function = @movmean;
	configuration.f0Window = 10;
	configuration.f1Function = @movstd;
	configuration.f1Window = 10;
	configuration.peaksLowpassFrequency = 0.2;
	configuration.thresholdingFunction = @std;
	configuration.thresholdFactor = 0.10;
	FPA(inputFile, configuration);
```

## Version History
* Initial Release: Library and example code

## License
© 2019 [Leonardo Molina][Leonardo Molina]

This project is licensed under the [GNU GPLv3 License][LICENSE.md].

[Leonardo Molina]: https://github.com/leomol
[MATLAB]: https://www.mathworks.com/downloads/
[LICENSE.md]: LICENSE.md