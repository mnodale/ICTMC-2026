# Code

This folder contains R code used to identify and review potential outliers in continuous clinical trial data.

## `detect_outliers.R`

Defines the `detect_outliers()` helper function used to flag potentially unusual continuous values.

The function expects data in long format, with:

- `vars`: variable name
- `value`: numeric value

For each variable, it calculates summary statistics and three standardized measures of dispersion:

- standard deviation-based z-score
- median absolute deviation (MAD) score
- interquartile range (IQR) score

Potential outliers can then be identified using any of these three approaches. The default method is MAD, providing a robust alternative to the usual mean/standard deviation approach.

The threshold is controlled by `k`; for example, `k = 0.995` corresponds to the two-sided standard-normal cut-off used in the accompanying QC code.

## `QC_outliers.R`

Applies the outlier-detection procedure to clinical trial data as part of data quality control.

The script:

- prepares continuous variables from several clinical data sources;
- examines their distributions;
- applies SD-, IQR-, and MAD-based outlier detection;
- uses the MAD-based results as the primary set of potential outliers for review;
- creates a catalogue of the variables and data points checked;
- exports potential outliers to an Excel workbook for review and data-query generation; and
- produces a text summary documenting the data reviewed and the number of potential outliers identified.

The code was developed as a pragmatic approach to identifying implausible or unusual values for clinical data quality control. Values flagged by the algorithm are intended for review rather than automatic exclusion or correction.

Marianna Nodale  
2026
