# ROCI Package

**Repository Type:** `roci` (R Package/Tool)

## Overview

ROCI is an R package that provides a set of functions to calculate and visualize the Receiver Operating Characteristic (ROC) curve, the Area Under the Curve (AUC), and other statistical metrics for binary classification models where the two classes are defined by an interval (i.e., two thresholds).

Unlike traditional ROC analysis that uses a single threshold, ROCI implements an algorithm for interval-based classification, where predictions are made based on whether values fall inside or outside a specified interval. This approach is particularly useful in medical and diagnostic applications where normal ranges are defined by lower and upper bounds.

## Key Features

- **Interval-based ROC Analysis**: Calculate ROC curves for binary classification using two cut-off values (intervals) instead of a single threshold
- **AUC Calculation**: Compute the Area Under the Curve for interval-based ROC curves
- **Visualization**: Generate publication-ready ROC curve plots with customizable features
- **Youden Index**: Calculate and highlight optimal intervals based on the Youden Index (J = sensitivity + specificity - 1)
- **Comprehensive Metrics**: Provides sensitivity, specificity, true positives, true negatives, false positives, and false negatives for all possible intervals

## Installation

You can install the development version of ROCI from GitHub:

```r
# Install devtools if not already installed
install.packages("devtools")

# Install ROCI package
devtools::install_github("Telehealth-Competence-Center-GmbH/ROCI_package")
```

## Usage

### Basic Example

```r
library(ROCI)

# Generate example data
set.seed(123)
ph_values <- rnorm(20, 7, 2)
med_class <- ifelse(ph_values >= quantile(ph_values, 0.25) & 
                    ph_values <= quantile(ph_values, 0.75),
                    sample(c(0, 1), prob = c(0.3, 0.7)),
                    sample(c(0, 1), prob = c(0.8, 0.2)))

# Calculate ROCI curve
roci_results <- roci(ph_values, med_class)

# Calculate AUC
auc_value <- auc_roci(roci_results)
print(paste("AUC:", round(auc_value, 3)))

# Plot the ROCI curve
plot_roci(roci_results)
```

### Main Functions

- `roci()`: Calculates sensitivity, specificity, and Youden Index for all possible intervals
- `auc_roci()`: Computes the Area Under the Curve for the ROCI curve
- `plot_roci()`: Creates a visualization of the ROCI curve with optional annotations

### Parameters

The `roci()` function supports classification where the event of interest can be either:
- **Inside the interval** (`event = "inside"`): Class 1 is predicted when values fall within the interval
- **Outside the interval** (`event = "outside"`): Class 1 is predicted when values fall outside the interval

## Dependencies

- dplyr
- ggplot2
- magrittr

## License

MIT License - see [LICENSE.md](LICENSE.md) for details

## Authors

- Dominik Thiele (dthiele@tcc-clinicalsolutions.de)
  - ORCID: 0000-0003-1073-0747

## References

For a full description of the ROCI algorithm, please refer to:
- Krannich 2023 (detailed reference to be added)

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Support

For questions or issues, please open an issue on the [GitHub repository](https://github.com/Telehealth-Competence-Center-GmbH/ROCI_package).
