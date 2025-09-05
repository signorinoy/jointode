
<!-- README.md is generated from README.Rmd. Please edit that file -->

# JointODE

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/ziyangg98/JointODE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ziyangg98/JointODE/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/ziyangg98/JointODE/graph/badge.svg)](https://app.codecov.io/gh/ziyangg98/JointODE)

<!-- badges: end -->

The **JointODE** package provides a unified framework for joint modeling
of longitudinal biomarker measurements and time-to-event outcomes using
ordinary differential equations (ODEs). This approach enables the
simultaneous analysis of biomarker trajectories and their impact on
survival outcomes.

## Model Setup

### Longitudinal Model

The observed biomarker measurements are modeled as:
$$V_{ij}=m_i(T_{ij})+b_i+\varepsilon_{ij},\quad i=1,\ldots,n,\quad j=1,\ldots,n_i$$

where:

- $V_{ij}$: Observed biomarker value for subject $i$ at time $T_{ij}$
- $m_i(t)$: True underlying biomarker trajectory
- $b_i\sim\mathcal{N}(0,\sigma_{b}^{2})$: Subject-specific random
  intercept
- $\varepsilon_{ij}\sim\mathcal{N}(0,\sigma_{e}^{2})$: Measurement error

The biomarker trajectory evolution is characterized by the following
second-order differential equation:

$$\ddot{m}_i(t) = f\big(m_i(t), \dot{m}_i(t), \mathbf{X}_i(t), t\big)$$

where
$f: \mathbb{R} \times \mathbb{R} \times \mathbb{R}^p \times \mathbb{R}^+ \to \mathbb{R}$
is a smooth function modeling the biomarker acceleration as a function
of its current value $m_i(t)$, velocity $\dot{m}_i(t)$, time-varying
covariates $\mathbf{X}_i(t) \in \mathbb{R}^p$, and time $t$.

### Survival Model

The hazard function incorporates biomarker dynamics:

$$\lambda_i(t) = \lambda_{0}(t)\exp\left[\mathbf{m}_i(t)^{\top}\boldsymbol{\alpha}+\mathbf{W}_i^{\top}\boldsymbol{\phi}+b_{i}\right]$$

where:

- $\lambda_{0}(t)$: Baseline hazard (e.g., Weibull, piecewise constant)
- $\mathbf{m}_i(t)=\left(m_i(t), \dot{m}_i(t)\right)^{\top}$: Biomarker
  value and velocity
- $\boldsymbol{\alpha}=(\alpha_1, \alpha_2)^{\top}$: Association
  parameters for value and velocity
- $\mathbf{W}_i$: Baseline covariates with coefficients
  $\boldsymbol{\phi}$
- $b_i$: Subject-specific random intercept

For detailed mathematical derivations including ODE formulation,
likelihood construction, and EM algorithm specifics, see the [technical
documentation](http://gongziyang.com/JointODE/articles/technical-details.html).

## Installation

You can install the development version of JointODE from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("ziyangg98/JointODE")
```

## Example

Here’s a basic example demonstrating typical usage:

``` r
library(JointODE)
#>
#> Attaching package: 'JointODE'
#> The following object is masked from 'package:stats':
#>
#>     simulate

# Load example dataset
data(sim)

# Fit joint ODE model
fit <- JointODE(
  longitudinal_formula = sim$formulas$longitudinal,
  longitudinal_data = sim$data$longitudinal_data,
  survival_formula = sim$formulas$survival,
  survival_data = sim$data$survival_data,
  parallel = TRUE
)

# Model summary
summary(fit)
#>
#> Call:
#> JointODE(longitudinal_formula = sim$formulas$longitudinal, longitudinal_data = sim$data$longitudinal_data,
#>     survival_formula = sim$formulas$survival, survival_data = sim$data$survival_data,
#>     parallel = TRUE)
#>
#> Variance components:
#> sigma_e sigma_b
#> 0.09686 0.09890
#>
#> Fixed effects:
#>                     Estimate Std. Error z value Pr(>|z|)
#> baseline:1         -3.183809   0.850407  -3.744 0.000181 ***
#> baseline:2         -3.056206   0.960422  -3.182 0.001462 **
#> baseline:3         -2.076161   0.698387  -2.973 0.002951 **
#> baseline:4         -1.482894   0.597729  -2.481 0.013106 *
#> baseline:5         -1.971558   0.619057  -3.185 0.001449 **
#> baseline:6         -1.770042   0.701663  -2.523 0.011648 *
#> baseline:7         -2.332177   1.456576  -1.601 0.109347
#> baseline:8         -0.387564   1.518789  -0.255 0.798584
#> baseline:9         -0.075105   1.393381  -0.054 0.957014
#> hazard:alpha1       0.593704   0.224000   2.650 0.008038 **
#> hazard:alpha2       0.919493   0.424881   2.164 0.030455 *
#> hazard:phi1         0.881763   0.149444   5.900 3.63e-09 ***
#> hazard:phi2        -1.471150   0.164842  -8.925  < 2e-16 ***
#> longitudinal:beta1 -1.013110   0.013813 -73.346  < 2e-16 ***
#> longitudinal:beta2 -0.604980   0.018811 -32.162  < 2e-16 ***
#> longitudinal:beta3 -0.001932   0.006591  -0.293 0.769465
#> longitudinal:beta4 -0.812085   0.012918 -62.866  < 2e-16 ***
#> longitudinal:beta5  0.500869   0.008569  58.449  < 2e-16 ***
#> longitudinal:beta6  0.405590   0.006437  63.014  < 2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#>
#> ---
#> Log-likelihood: 918.1049   AIC: -1794.21   BIC: -1739.501
#> N = 100  Convergence: EM algorithm converged after 21 iterations

# Generate predictions
predictions <- predict(fit, times = seq(0, 10, by = 0.25))
```

## Visualization

    #>
    #> Attaching package: 'dplyr'
    #> The following objects are masked from 'package:stats':
    #>
    #>     filter, lag
    #> The following objects are masked from 'package:base':
    #>
    #>     intersect, setdiff, setequal, union

<img src="man/figures/README-visualization-1.png" width="100%" />

## Code of Conduct

Please note that the JointODE project is released with a [Contributor
Code of Conduct](http://gongziyang.com/JointODE/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
