## Requirements

- The package requires the installation of the JAGS software. Click to [download JAGS](https://sourceforge.net/projects/mcmc-jags/).

## New Features

- **Multi-core Processing**: The package now supports parallel processing for improved performance. MCMC chains and cross-validation folds can run in parallel across multiple CPU cores.
- **Cross-platform Compatibility**: Optimized for Unix/Linux, macOS, and Windows systems.
- **Backward Compatibility**: All existing code continues to work without modification.

## Installation

- This package is not currently on cran so you can download from Github. Make sure to have the `devtools` package installed and then execute the following:

```
devtools::install_github("ncahill89/BTFr")
```

You can then load the BTF package using the `library` function.

```{r}
library(BTFr)
```

## Getting started

See [Vignette](https://github.com/ncahill89/vignettes/blob/master/BTF.md).

### Parallel Processing

The package now supports multi-core processing for improved performance:

```r
# Default behavior uses all available cores minus 1
modern_mod <- run_modern(modern_elevation = NJ_modern_elevation,
                         modern_species = NJ_modern_species)

# Specify number of cores
modern_mod <- run_modern(modern_elevation = NJ_modern_elevation,
                         modern_species = NJ_modern_species,
                         n_cores = 4)

# Disable parallel processing for original behavior
modern_mod <- run_modern(modern_elevation = NJ_modern_elevation,
                         modern_species = NJ_modern_species,
                         parallel = FALSE)
```

For detailed information about parallel processing features, see `PARALLEL_OPTIMIZATION.md`.

## Workflow Example (parallel processing)

see [test.R](test.R) for a complete workflow example.

### Performance Benchmarks

**Tested on 16-core laptop:**

- **Modern calibration**: ~5.5 minutes
- **Core reconstruction**: ~8.5 minutes
- **10-fold cross-validation**: ~14 minutes
- **Total workflow**: ~28 minutes

### Memory and Resource Tips

```r
# For systems with limited memory, reduce core usage:
run_modern(parallel = TRUE, n_cores = 2)

# For very large datasets, consider sequential processing:
run_modern(parallel = FALSE)

# Monitor system resources:
cat("Available cores:", detectCores(), "\n")
cat("Memory usage:", round(object.size(modern_mod) / 1024^2, 1), "MB\n")
```

## Data

There are options to use the package default which contain data for New Jersey, USA. Alternatively, you can supply your own data. When supplying data, use the package defaults as templates for formatting. You'll find more details in [the vignette](https://github.com/ncahill89/vignettes/blob/master/BTF.md).
