
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `R`/`unihtee`

> Univariate Heterogeneous Treatment Effect Estimation

**Author:** [Philippe Boileau](https://pboileau.ca/)

<!-- badges: start -->

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
<!-- badges: end -->

------------------------------------------------------------------------

`unicate` provides tools for uncovering treatment effect modifiers in
high-dimensional data. Treatment effect modification is defined using
variable importance parameters based on absolute and relative effects.
The variable importance parameters are model-agnostic, meaning they do
not depend on specific regression pocedures (for example, the variable
importance measures of Random Forests). Inference is performed about
these variable importance measures using nonparametric estimators. Users
may use one-step or targeted maximum likelihood estimators. Under
general conditions, these estimators are unbiased and efficient.

## Installation

The package may be installed from GitHub using
[`remotes`](https://CRAN.R-project.org/package=remotes):

``` r
remotes::install_github("insightsengineering/unihtee")
```

`unihtee` is under active development. Check back often for updates.

## Usage

`unihtee()` is the only user-facing function. It can be used to perform
inference about the treatment effect modification variable importance
parameters. These parameters are defined for use with continuous, binary
and time-to-event outcomes with binary exposure variables. Variable
importance parameters are defined on both the absolute and relative
effect scales.

## Example

We simulate some observational study data that contains ten confounders,
of which are two treatment effect modifiers. We then perform inference
about the absolute treatment effect modifier variable importance
parameter.

``` r
library(unihtee)
library(MASS)
library(data.table)
library(sl3)

set.seed(510)

## create the dataset
n_obs <- 500
w <- mvrnorm(n = n_obs, mu = rep(0, 10), Sigma = diag(10))
confounder_names <- paste0("w_", seq_len(10))
colnames(w) <- confounder_names
a <- rbinom(n = n_obs, size = 1, prob = plogis(w[, 1] + w[, 2]))
y <- rnorm(n = n_obs, mean = w[, 1] + w[, 2] + a * w[, 3] - a * w[, 4])
dt <- as.data.table(cbind(w, a, y))

## targeted maximum likelihood estimates and testing procedure
unihtee(
  data = dt,
  confounders = confounder_names,
  modifiers = confounder_names,
  exposure = "a",
  outcome = "y",
  outcome_type = "continuous",
  effect = "absolute",
  estimator = "tmle"
)
#>     modifier     estimate        se           z      p_value    ci_lower
#>  1:      w_3  1.044592804 0.1613285  6.47494319 9.484769e-11  0.72838896
#>  2:      w_4 -0.869002514 0.1492388 -5.82289742 5.783606e-09 -1.16151066
#>  3:      w_8  0.137803254 0.1137965  1.21096238 2.259098e-01 -0.08523784
#>  4:      w_1  0.115258422 0.1160997  0.99275414 3.208298e-01 -0.11229692
#>  5:      w_9  0.124150185 0.1300374  0.95472664 3.397160e-01 -0.13072315
#>  6:     w_10 -0.097928234 0.1356976 -0.72166517 4.705004e-01 -0.36389554
#>  7:      w_6  0.054845105 0.1159964  0.47281713 6.363437e-01 -0.17250792
#>  8:      w_2 -0.064478504 0.1767632 -0.36477331 7.152806e-01 -0.41093441
#>  9:      w_7 -0.014704981 0.1485331 -0.09900136 9.211372e-01 -0.30582989
#> 10:      w_5  0.001500152 0.1103752  0.01359138 9.891560e-01 -0.21483526
#>       ci_upper  p_value_fdr
#>  1:  1.3607966 9.484769e-10
#>  2: -0.5764944 2.891803e-08
#>  3:  0.3608444 6.794319e-01
#>  4:  0.3428138 6.794319e-01
#>  5:  0.3790235 6.794319e-01
#>  6:  0.1680391 7.841673e-01
#>  7:  0.2821981 8.941008e-01
#>  8:  0.2819774 8.941008e-01
#>  9:  0.2764199 9.891560e-01
#> 10:  0.2178356 9.891560e-01
```

## Issues

If you encounter any bugs or have any specific feature requests, please
[file an issue](https://github.com/insightsengineering/unihtee/issues).

## Contributions

Contributions are very welcome. Interested contributors should consult
our [contribution
guidelines](https://github.com/insightsengineering/unihtee/blob/master/.github/CONTRIBUTING.md)
prior to submitting a pull request.

## License

The contents of this repository are distributed under the Apache 2.0
license. See the
[`LICENSE.md`](https://github.com/insightsengineering/unihtee/blob/main/LICENSE.md)
and
[`LICENSE`](https://github.com/insightsengineering/unihtee/blob/main/LICENSE)
files for details.
