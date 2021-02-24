## Installation 
Download the BTF package and use `devtools::install("BTF")` to install the package locally on your machine. Make sure that your working directory is set to where you have saved the BTF package. 

Note the package will not work unless you have the JAGS (Just Another Gibbs Sampler) software installed on your machine. You can download and install JAGS from [here](https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Windows/) for Windows and [here](https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Mac%20OS%20X/) for MAC. 

```{r}
devtools::install("BTF")
```

You can then load the BTF package using the `library` function. 
```{r}
library(BTF)
```

## About
This package can be used to run the Bayesian Transfer Function (BTF). 
There are options to use the package default which contains data for New Jersey, USA. 
Alternatively, you can supply your own data. When supplying data, use the package defaults as templates for formatting. You'll find the vignette [here](https://rpubs.com/ncahill_stat/720476).

