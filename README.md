
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
  survival_formula = sim$formulas$survival,
  longitudinal_data = sim$data$longitudinal_data,
  survival_data = sim$data$survival_data
)

# Model summary
summary(fit)
#>
#> Call:
#> JointODE(longitudinal_formula = sim$formulas$longitudinal, survival_formula = sim$formulas$survival,
#>     longitudinal_data = sim$data$longitudinal_data, survival_data = sim$data$survival_data)
#>
#> Data Descriptives:
#> Longitudinal Process            Survival Process
#> Number of Observations: 1529    Number of Events: 73 (73%)
#> Number of Subjects: 100
#>
#>        AIC        BIC     logLik
#>  -2592.179  -2545.286   1314.090
#>
#> Coefficients:
#> Longitudinal Process: Second-Order ODE Model
#>              Estimate Std. Error z value Pr(>|z|)
#> v_value     -0.600890   0.003737 -160.80   <2e-16 ***
#> v_slope     -0.400776   0.006375  -62.87   <2e-16 ***
#> (Intercept)  0.294265   0.002673  110.07   <2e-16 ***
#> x1          -0.800519   0.005831 -137.28   <2e-16 ***
#> x2           0.497675   0.003858  129.00   <2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#>
#> Survival Process: Proportional Hazards Model
#>         Estimate Std. Error z value Pr(>|z|)
#> v_value   0.4462     0.1233   3.617 0.000298 ***
#> v_slope   0.9769     0.3416   2.860 0.004243 **
#> w1        0.5006     0.1364   3.671 0.000241 ***
#> w2       -0.7852     0.1241  -6.329 2.46e-10 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#>
#> Baseline Hazard: B-spline with 7 basis functions
#> (Coefficients range: [-3.646, -0.640] )
#>
#> Variance Components:
#>               StdDev
#> Random Effect       0.101471
#> Residual            0.098552
#>
#> Model Diagnostics:
#> C-index (Concordance): 0.584
#> Convergence: EM algorithm converged after 32 iterations

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
