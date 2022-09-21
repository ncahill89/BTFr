## Requirements

  - The package requires the installation of the JAGS software. Click to [download JAGS](https://sourceforge.net/projects/mcmc-jags/).

## Installation

  - This package is not currently on cran so you can download from Github. Make sure to have the `devtools` package installed and then execute the following: 

```
devtools::install_github("ncahill89/BTF")
```

You can then load the BTF package using the `library` function. 

```{r}
library(BTF)
```

## Getting started

See [Vignette](https://rpubs.com/ncahill_stat/720476).

## Data

There are options to use the package default which contain data for New Jersey, USA. Alternatively, you can supply your own data. When supplying data, use the package defaults as templates for formatting. You'll find more details in [the vignette](https://github.com/ncahill89/vignettes/blob/master/BTF.md).

