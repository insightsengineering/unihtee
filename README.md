
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `R`/`unihtee`

> Univariate Heterogeneous Treatment Effect Estimation

**Author:** [Philippe Boileau](https://pboileau.ca/)

<!-- badges: start -->
<!-- badges: end -->

------------------------------------------------------------------------

`unicate` provides tools for uncovering treatment effect modifiers in
high-dimensional data. Treatment effect modification is defined using
variable importance parameters based on risk differences or relative
risks. The variable importance parameters are model-agnostic, meaning
they do not depend on specific regression pocedures (for example, the
variable importance measures of Random Forests). Inference is performed
about these variable importance measures using nonparametric estimators.
Users may use one-step or targeted maximum likelihood estimators. Under
general conditions, these estimators are double-robust, linear and
efficient.

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
importance parameters are defined on both the risk difference and
relative risk scales.

## Example

We simulate some observational data that contains ten confounders, of
which are two treatment effect modifiers. We then perform inference
about the risk-difference-based variable importance parameter.

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
  risk_type = "risk difference",
  estimator = "tmle"
)
#>     modifier     estimate        se           z      p_value    ci_lower
#>  1:      w_3  1.108005609 0.1613233  6.86823260 6.500134e-12  0.79181204
#>  2:      w_4 -0.839689658 0.1492369 -5.62655519 1.838440e-08 -1.13219398
#>  3:      w_8  0.188515268 0.1137877  1.65672738 9.757461e-02 -0.03450871
#>  4:      w_1  0.177774878 0.1160698  1.53162049 1.256161e-01 -0.04972192
#>  5:     w_10 -0.109936815 0.1356947 -0.81017792 4.178379e-01 -0.37589834
#>  6:      w_9  0.103756298 0.1300309  0.79793593 4.249077e-01 -0.15110420
#>  7:      w_6  0.076346514 0.1159920  0.65820496 5.104064e-01 -0.15099782
#>  8:      w_2 -0.082664876 0.1767372 -0.46772776 6.399793e-01 -0.42906971
#>  9:      w_7 -0.018588007 0.1485309 -0.12514570 9.004082e-01 -0.30970862
#> 10:      w_5  0.001748801 0.1103689  0.01584504 9.873580e-01 -0.21457433
#>       ci_upper  p_value_fdr
#>  1:  1.4241992 6.500134e-11
#>  2: -0.5471853 9.192199e-08
#>  3:  0.4115392 3.140403e-01
#>  4:  0.4052717 3.140403e-01
#>  5:  0.1560247 7.081795e-01
#>  6:  0.3586168 7.081795e-01
#>  7:  0.3036908 7.291521e-01
#>  8:  0.2637400 7.999741e-01
#>  9:  0.2725326 9.873580e-01
#> 10:  0.2180719 9.873580e-01
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
