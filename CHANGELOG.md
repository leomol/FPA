## Changelog
* 2022-12-16
	- Added option for detecting peaks using prominence or height (default). When prominence is selected, threshold is not plotted.
	- Changed the way thresholds are provided; using a raw value or a function that operates on input data.
	- Changed default peak separation from 2.0 to 0.5.
	- Removed rendering of negative peaks.
	- Removed examples. Will re-add on a request basis.
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