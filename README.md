
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R/`unihtee`

> Univariate Heterogeneous Treatment Effect Estimation

**Author:** [Philippe Boileau](https://pboileau.ca/)

<!-- badges: start -->

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
<!-- badges: end -->

------------------------------------------------------------------------

`unihtee` provides tools for uncovering treatment effect modifiers in
high-dimensional data. Treatment effect modification is defined using
variable importance parameters based on absolute and relative effects.
Inference is performed about these variable importance measures using
nonparametric estimators. Users may use one-step or targeted maximum
likelihood estimators. Under general conditions, these estimators are
unbiased and efficient.

Additional details about this methodology are provided in Boileau et al.
(2022), Boileau et al. (2025), and in the package’s
[vignette](https://insightsengineering.github.io/unihtee/main/articles/using-unihtee.html).

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
parameters. These parameters are defined for data-generating processes
with continuous, binary and time-to-event outcomes with binary exposure
variables. Variable importance parameters based on absolute and relative
effects are available. Details are provided in the vignette.

## Example

We simulate some observational study data that contains ten
pre-treatment covariates, of which are two treatment effect modifiers.
We then perform inference about the absolute treatment effect modifier
variable importance parameter, which is inspired by the average
treatment effect.

``` r
library(unihtee)
library(MASS)
library(data.table)
library(sl3)

set.seed(510)

# create the dataset
n_obs <- 500
w <- mvrnorm(n = n_obs, mu = rep(0, 10), Sigma = diag(10))
confounder_names <- paste0("w_", seq_len(10))
colnames(w) <- confounder_names
a <- rbinom(n = n_obs, size = 1, prob = plogis(w[, 1] + w[, 2]))
y <- rnorm(n = n_obs, mean = w[, 1] + w[, 2] + a * w[, 3] - a * w[, 4])
dt <- as.data.table(cbind(w, a, y))

# targeted maximum likelihood estimates and testing procedure
unihtee_output <- unihtee(
  data = dt,
  confounders = confounder_names,
  modifiers = confounder_names,
  exposure = "a",
  outcome = "y",
  outcome_type = "continuous",
  effect = "absolute",
  estimator = "tmle"
)

# plot table of estimated TEM-VIPs
unihtee_output$temvip_inference_tbl
#>     modifier     estimate        se           z      p_value    ci_lower
#>       <fctr>        <num>     <num>       <num>        <num>       <num>
#>  1:      w_3  1.044592804 0.1599527  6.53063474 6.549161e-11  0.73108547
#>  2:      w_4 -0.869002514 0.1505155 -5.77351005 7.763697e-09 -1.16401281
#>  3:      w_8  0.137803254 0.1138565  1.21032372 2.261547e-01 -0.08535554
#>  4:      w_1  0.115258422 0.1168820  0.98610900 3.240796e-01 -0.11383036
#>  5:      w_9  0.124150185 0.1295567  0.95826884 3.379272e-01 -0.12978102
#>  6:     w_10 -0.097928234 0.1356345 -0.72200119 4.702937e-01 -0.36377176
#>  7:      w_6  0.054845105 0.1157812  0.47369617 6.357166e-01 -0.17208602
#>  8:      w_2 -0.064478504 0.1761998 -0.36593963 7.144101e-01 -0.40983019
#>  9:      w_7 -0.014704981 0.1485610 -0.09898279 9.211519e-01 -0.30588450
#> 10:      w_5  0.001500152 0.1100365  0.01363321 9.891226e-01 -0.21417147
#>       ci_upper  p_value_fdr
#>          <num>        <num>
#>  1:  1.3581001 6.549161e-10
#>  2: -0.5739922 3.881849e-08
#>  3:  0.3609620 6.758544e-01
#>  4:  0.3443472 6.758544e-01
#>  5:  0.3780814 6.758544e-01
#>  6:  0.1679153 7.838229e-01
#>  7:  0.2817762 8.930127e-01
#>  8:  0.2808732 8.930127e-01
#>  9:  0.2764745 9.891226e-01
#> 10:  0.2171718 9.891226e-01
```

## Issues

If you encounter any bugs or have any specific feature requests, please
[file an issue](https://github.com/insightsengineering/unihtee/issues).

## Contributions

Contributions are very welcome. Interested contributors should consult
our [contribution
guidelines](https://github.com/insightsengineering/unihtee/blob/master/.github/CONTRIBUTING.md)
prior to submitting a pull request.

## Citation

To cite `unihtee` and the papers introducing the underlying framework,
use the following BibTeX entries:

    @manual{unihtee,
      title = {unihtee: Univariate Heterogeneous Treatment Effect Estimation},
      author = {Philippe Boileau},
      note = {R package version 0.0.1}
    }

    @article{boileau2025,
      title = {A Nonparametric Framework for Treatment Effect Modifier Discovery in High Dimensions},
      author = {Boileau, Philippe and Leng, Ning and Hejazi, Nima S and {van der Laan}, Mark and Dudoit, Sandrine},
      year = {2025},
      journal = {Journal of the Royal Statistical Society Series B: Statistical Methodology},
      volume = {87},
      number = {1},
      pages = {157--185},
      issn = {1369-7412},
      doi = {10.1093/jrsssb/qkae084}
    }

    @article{boileau2022,
      author = {Boileau, Philippe and Qi, Nina Ting and van der Laan, Mark J and Dudoit, Sandrine and Leng, Ning},
      title = {A flexible approach for predictive biomarker discovery},
      journal = {Biostatistics},
      year = {2022},
      month = {07},
      issn = {1465-4644},
      doi = {10.1093/biostatistics/kxac029},
      url = {https://doi.org/10.1093/biostatistics/kxac029}
    }

## License

The contents of this repository are distributed under the Apache 2.0
license. See the
[`LICENSE.md`](https://github.com/insightsengineering/unihtee/blob/main/LICENSE.md)
and
[`LICENSE`](https://github.com/insightsengineering/unihtee/blob/main/LICENSE)
files for details.

## References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-boileau2023" class="csl-entry">

Boileau, Philippe, Ning Leng, Nima S Hejazi, Mark van der Laan, and
Sandrine Dudoit. 2025. “A Nonparametric Framework for Treatment Effect
Modifier Discovery in High Dimensions.” *Journal of the Royal
Statistical Society Series B: Statistical Methodology*.
<https://doi.org/10.1093/jrsssb/qkae084>.

</div>

<div id="ref-boileau2022" class="csl-entry">

Boileau, Philippe, Nina Ting Qi, Mark J van der Laan, Sandrine Dudoit,
and Ning Leng. 2022. “<span class="nocase">A flexible approach for
predictive biomarker discovery</span>.” *Biostatistics*, July.
<https://doi.org/10.1093/biostatistics/kxac029>.

</div>

</div>
