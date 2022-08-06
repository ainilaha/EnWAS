# Environment‐Wide Association Study (EnWAS)

EnWAS is .................................

## Installation

You can install the development version of BioImpute from [GitHub](https://github.com/) with:

``` {r}
# install.packages("devtools")
devtools::install_github("ccb-hms/EnWAS")
```

## Example

This is a basic example which shows you how to solve a common problem:

```{r}
library(EnWAS)
data(exposure_vars)
data(nhanes)
linear_model <- 'BMXWAIST ~ RIDAGEYR*RIAGENDR + BMXBMI'
linear_res <- enwas(linear_model, exposure_vars, nhanes)
```
