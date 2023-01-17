
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `R`/`unihtee`

> Univariate Heterogeneous Treatment Effect Estimation

**Author:** [Philippe Boileau](https://pboileau.ca/)

<!-- badges: start -->
<!-- badges: end -->

------------------------------------------------------------------------

`unicate` provides tools for uncovering treatment effect modifiers in
high-dimensional data. Treatment effect modification is defined using
variable importance parameters based on absolute and relative scales.
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
scales.

## Example

We simulate some observational study data that contains ten confounders,
of which are two treatment effect modifiers. We then perform inference
about the absolute scale variable importance parameter.

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
  scale = "absolute",
  estimator = "tmle"
)
#>     modifier    estimate        se          z      p_value    ci_lower
#>  1:      w_3  1.19657424 0.1613233  7.4172461 1.196820e-13  0.88038067
#>  2:      w_4 -0.89218696 0.1492369 -5.9783268 2.254412e-09 -1.18469128
#>  3:      w_8  0.19677508 0.1137877  1.7293170 8.375237e-02 -0.02624889
#>  4:      w_1  0.19352745 0.1160698  1.6673368 9.544745e-02 -0.03396935
#>  5:     w_10 -0.11398276 0.1356947 -0.8399945 4.009115e-01 -0.37994429
#>  6:      w_9  0.10749149 0.1300309  0.8266613 4.084291e-01 -0.14736901
#>  7:      w_6  0.07756350 0.1159920  0.6686969 5.036888e-01 -0.14978083
#>  8:      w_2 -0.08609331 0.1767372 -0.4871262 6.261689e-01 -0.43249814
#>  9:      w_7 -0.01905859 0.1485309 -0.1283140 8.979005e-01 -0.31017921
#> 10:      w_5  0.00176980 0.1103689  0.0160353 9.872062e-01 -0.21455333
#>       ci_upper  p_value_fdr
#>  1:  1.5127678 1.196820e-12
#>  2: -0.5996826 1.127206e-08
#>  3:  0.4197991 2.386186e-01
#>  4:  0.4210242 2.386186e-01
#>  5:  0.1519788 6.807151e-01
#>  6:  0.3623520 6.807151e-01
#>  7:  0.3049078 7.195555e-01
#>  8:  0.2603115 7.827111e-01
#>  9:  0.2720620 9.872062e-01
#> 10:  0.2180929 9.872062e-01
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
