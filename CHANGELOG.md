## Changelog
* v2.0.0
	* 2023-11-10
		- Changed version from `1.0.0` to `2.0.1` in the citation file.
		- Changed `plotTraces` to `plotTrace` in the documentation.
	* 2023-11-09
		- New changes saved to a new branch.
		- Packaged data loaders.
		- Renamed data properties to match the name of the corresponding function that generates it.
		- Analysis steps are part of the configuration.
		- Outcomes of the analysis are read-only.
		- Added standardization step in the analysis and plots.
		- Added 95% confidence interval to PowerSpectrumPlot and uses the Welch's method.
		- Added helper functions to customize analysis steps.
		- Window size and event times moved from the configuration to an argument of the plotting or exporting functions.
		- Exporting data may now occur in rows or columns.
		- Examples are described in the online documentation.
		- Removed deprecated data and example files.
* v1.0.0
	* 2023-10-12
		- Configuration is now reflected as members of the FPA object.
		- Changed the signature of multiple parameters to provide flexibility with the analysis steps.
			- The following parameters are affected: `configuration.threshold`, `configuration.signalCorrection` and `configuration.referenceCorrection`. Now these functions expect `dff` as argument, `configuration.threshold` expect a single number as an output argument, whereas `signal` and `referenceCorrection` expect an array of the same size of `fpa.dff` as output.
		- User provided epochs are now used _as is_ without consideration of `edgeId` or `artifactFreeId`.
		- Removed `thresholdEpochs` as the user may incorporate this information in the aformetioned function.
		- Added `signalArtifactFree` and `referenceArtifactFree` to `fpa`'s properties.
	* 2023-04-12
		- Added option to modify the default baseline correction operation (signal bleaching is subtracted from signal and reference bleaching is subtracted from reference).
		- Changed the threshold function signature; arguments to calling function are now fpa and data.
	* 2023-04-05
		- Added option to provide separate baseline epochs for signal and reference.
	* 2023-02-03
		- Added h5 loader, which works for Doric lock-in datasets by default.
		- Fix an out of bounds bug that could happen for peaks at the beginning or end of a trace.
		- Removing time points later than the original time series.
	* 2022-12-16
		- Added option for detecting peaks using prominence or height (default). When prominence is selected, threshold is not plotted.
		- Changed the way thresholds are provided; using a raw value or a function that operates on input data.
		- Changed default peak separation from 2.0 to 0.5.
		- Removed rendering of negative peaks.
		- Removed examples. Will re-add on a request basis.
		- Event and peak windows can be defined as a range or as a vector.
		- Added option to plot event trigger average.
		- Added option normalize event data prior to plotting (undocumented).
		- Added a reader for Deep Lab Cut generated files (undocumented).
		- Added normalization options for peaks.
		- Event and peak are normalized to values left of the peak by default.
	* 2022-08-31
		- Updated `exampleRA.m` to show how to normalize.
		- Documented code and `readme.md`
	* 2022-06-30
		- Added setup.
		- Added exponential decay fit to trigger plots.
		- Removed smooth dff plots and object fields.
		- Thresholding value is calculated on both artifact and edge free data.
		- Removed `peaksLowpassFrequency`:
			- Peak detection is done directly in dff with minimum distance given by peakSeparation (convert with `1.0 / peaksLowpassFrequency`).
		- Filtering includes data edges.
		- Fitting reference to signal uses all data.
		- Added event window for event triggers.
		- Added normalization options for events.
		- Refactored triggered window to peak window which can now be a range.
		- Removed start and stop triggers; using events instead.
		- Fixed epoch bands not drawing epochs with Inf values.
	* 2022-06-08
		- Added custom event triggers (similar to start, stop, and peak triggers).
		- Added custom epochs to calculate peak threshold.
		- All triggered data exported includes the column ConditionName.
	* 2022-05-10
		- Export f.
		- Added batch processing scripts to traverse project folders.
	* 2021-11-17
		- Patched `loadData.m` to remove columns with all NaN values and rows with any NaN values.
		- Added column `normalizedArea` (redundantly) to exported AUC file.
	* 2021-11-01
		- Process a project folder as an example.
		- First release.
	* 2021-03-25
		- Added baseline plot of reference for figure 1.1.
		- Figure 1.1 shows resampled raw data (before was corrected data).
	* 2021-03-15
		- Added heatmaps of raw triggers.
		- Plots are now methods.
	* 2021-03-02
		- Added thresholding options.
		- Default thresholding changed from k * mad(x) + mean(x) to k * mad(x) + median(x).
	* 2021-03-01
		- FPA changed from function to object.
		- Added export method.
		- Added epoch-start and epoch-stop to trigger averages.
		- Added BORIS reader.
	* 2021-02-25
		- Area under the curve is normalized and calculated with trapz.
		- Default normalization changed to "modified z-score".
		- Significant changes to normalization API.
	* 2021-01-15
		- Reference is fitted to signal.
		- Band-pass replaced by low-pass.
		- GUI example supports artifact removal.
	* 2019-08-22
		- Library and example code.
		- Added loaders for multiple data acquisition systems.
	* 2019-02-20
		- Initial release.

## License
© 2018 [Leonardo Molina][Leonardo Molina]

This project is licensed under the [GNU GPLv3 License][LICENSE.md].

[Leonardo Molina]: https://github.com/leomol
[LICENSE.md]: LICENSE.md