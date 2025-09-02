
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
- $\mathbf{m}_i(t)=\left(m_i(t), \dot{m}_i(t), \ddot{m}_i(t)\right)^{\top}$:
  Biomarker value and derivatives
- $\boldsymbol{\alpha}=(\alpha_0, \alpha_1, \alpha_2)^{\top}$:
  Association parameters for value, velocity, and acceleration
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
  init = sim$parameters,
  parallel = TRUE
)

# Model summary
summary(fit)
#>
#> Call:
#> JointODE(longitudinal_formula = sim$formulas$longitudinal, longitudinal_data = sim$data$longitudinal_data,
#>     survival_formula = sim$formulas$survival, survival_data = sim$data$survival_data,
#>     init = sim$parameters, parallel = TRUE)
#>
#> Variance components:
#> sigma_e sigma_b
#> 0.10092 0.09506
#>
#> Fixed effects:
#>                     Estimate Std. Error z value Pr(>|z|)
#> baseline:1         -3.867098   0.844898  -4.577 4.72e-06 ***
#> baseline:2         -2.554384   0.746438  -3.422 0.000621 ***
#> baseline:3         -2.469201   0.619648  -3.985 6.75e-05 ***
#> baseline:4         -2.070030   0.441657  -4.687 2.77e-06 ***
#> baseline:5         -1.974180   0.409197  -4.825 1.40e-06 ***
#> baseline:6         -1.816067   0.470518  -3.860 0.000114 ***
#> baseline:7         -1.657955   0.900964  -1.840 0.065739 .
#> baseline:8         -1.540116   1.267635  -1.215 0.224384
#> baseline:9         -1.487824   1.350214  -1.102 0.270497
#> hazard:alpha0       0.468010   0.068575   6.825 8.80e-12 ***
#> hazard:alpha1       0.006911   0.254299   0.027 0.978318
#> hazard:alpha2       0.301793   0.307082   0.983 0.325717
#> hazard:phi1         0.373701   0.083733   4.463 8.08e-06 ***
#> hazard:phi2        -0.640911   0.082446  -7.774 7.62e-15 ***
#> longitudinal:beta1 -0.989857   0.013367 -74.052  < 2e-16 ***
#> longitudinal:beta2 -1.976997   0.037872 -52.202  < 2e-16 ***
#> longitudinal:beta3  0.005200   0.006199   0.839 0.401579
#> longitudinal:beta4  1.980852   0.029870  66.317  < 2e-16 ***
#> longitudinal:beta5  0.990293   0.015023  65.920  < 2e-16 ***
#> longitudinal:beta6  0.493552   0.008160  60.481  < 2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#>
#> ---
#> Log-likelihood: 2042.664   AIC: -4041.328   BIC: -3968.765
#> N = 200  Convergence: EM algorithm converged after 6 iterations

# Generate predictions
predictions <- predict(fit, times = seq(0, 10, by = 0.25))
```

## Visualization

``` r
library(ggplot2)
library(dplyr)
#>
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#>
#>     filter, lag
#> The following objects are masked from 'package:base':
#>
#>     intersect, setdiff, setequal, union
library(tidyr)

# Prepare data
ids <- c(1, 5, 10, 15)
df <- lapply(ids, \(i) {
  pred <- predictions[[i]]
  obs <- filter(sim$data$longitudinal_data, id == i)
  data.frame(
    id = i, time = pred$times,
    biomarker_est = pred$biomarker, velocity_est = pred$velocity,
    acceleration_est = pred$acceleration, survival = pred$survival
  ) %>%
    left_join(select(obs, time, v,
      biomarker_true = biomarker,
      velocity_true = velocity, acceleration_true = acceleration
    ), by = "time")
}) %>% bind_rows()

# Biomarker plot
ggplot(df, aes(x = time)) +
  geom_line(aes(y = biomarker_est, color = "Estimated")) +
  geom_line(
    aes(y = biomarker_true, color = "True"), linetype = 2, na.rm = TRUE
  ) +
  geom_point(aes(y = v, color = "Observed"), alpha = 0.3, na.rm = TRUE) +
  facet_wrap(~id) +
  theme_minimal() +
  labs(y = "Biomarker", color = "")
```

<img src="man/figures/README-visualization-1.png" width="100%" />

``` r

# Dynamics plot
df %>%
  pivot_longer(matches("(velocity|acceleration)_(est|true)"),
    names_to = c("type", "source"), names_sep = "_"
  ) %>%
  filter(!is.na(value)) %>%
  ggplot(aes(x = time, y = value, color = type, linetype = source)) +
  geom_line() +
  facet_wrap(~id, scales = "free_y") +
  theme_minimal() +
  labs(y = "Value", color = "Type", linetype = "Source")
```

<img src="man/figures/README-visualization-2.png" width="100%" />

``` r

# Survival plot
ggplot(df, aes(x = time, y = survival, color = factor(id))) +
  geom_line(linewidth = 1.2) +
  theme_minimal() +
  ylim(0, 1) +
  labs(y = "Survival Probability", color = "Subject")
```

<img src="man/figures/README-visualization-3.png" width="100%" />

## Code of Conduct

Please note that the JointODE project is released with a [Contributor
Code of Conduct](http://gongziyang.com/JointODE/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
