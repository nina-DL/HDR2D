# HDR2D: Alternative Approaches for Estimating Highest-Density Regions

The R package **HDR2D** provides a framework for estimating highest-density regions in two dimensions. It uses different measures, including kernel density estimation and some generalized approaches based on neaighborhood measures. 

The package implements the measures described in the following paper: [Nina Deliu and Brunero Liseo. "Alternative Approaches for Estimating Highest-Density Regions". arXiv preprint arXiv:2401.00245 (2023)](https://arxiv.org/abs/2401.00245)


## Installation

You can install the *development* version of **HDR2D** from Github using the commands:

``` r
install.packages("devtools")
devtools::install_github("nina-DL/HDR2D")
```

## Example

This is a basic example that shows you the use the main function **HDR.2d()**.

```{r example}
library(HDR2D)

# Generate some bivariate data
R = matrix(c(1, 0.5, 0.5, 1), nrow = 2, ncol = 2, byrow = TRUE)
draws_2d = MASS::mvrnorm(n = 2000, mu = c(0, 0), Sigma = R)

# Estimate the HDR with KDE
HDR.2d(sn = draws_2d, measure = "KDE", coverage_prob = 0.95, build_plot = T)

# Estimate the HDR with KNN
HDR.2d(sn = draws_2d, measure = "KNN", coverage_prob = 0.95, build_plot = T)

# Estimate the HDR with Parametric CDF Copula
HDR.2d(sn = draws_2d, measure = "CDF.PCop", coverage_prob = 0.95, build_plot = T, margin_family = c('norm', 't'))

```
