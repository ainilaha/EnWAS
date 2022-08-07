# Environment‐Wide Association Study (EnWAS)

EnWAS is .................................

## Installation

You can install the development version of BioImpute from [GitHub](https://github.com/) with:

``` {r}
# install.packages("devtools")
devtools::install_github("ccb-hms/EnWAS")
```

## Examples

This is a basic example which shows you how to solve a common problem:

```{r}
library(EnWAS)
data(exposure_vars)
data(nhanes)
linear_model <- 'BMXWAIST ~ RIDAGEYR*RIAGENDR + BMXBMI'
linear_res <- enwas(linear_model, exposure_vars, nhanes)
```


Plot result

```{r}
forest_plot(linear_res$enwas_res)
```

Plot results

```{r}
ns_model <-
  'BMXWAIST ~ splines::ns(RIDAGEYR, knots = seq(30, 80, by = 10)) * RIAGENDR + splines::ns(BMXBMI,knots = seq(30, 80, by = 10))'
ns_res <- enwas(ns_model, exposure_vars, nhanes)

forest_plot_mult(
  list(
    linear = linear_res$enwas_res,
    ns = ns_res$enwas_res
  )
)
```





