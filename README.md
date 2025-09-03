
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
#> 0.09592 0.10462 
#> 
#> Fixed effects:
#>                     Estimate Std. Error z value Pr(>|z|)    
#> baseline:1         -2.084879   0.250521  -8.322  < 2e-16 ***
#> baseline:2         -2.070379   0.838190  -2.470  0.01351 *  
#> hazard:alpha0       0.675740   0.221085   3.056  0.00224 ** 
#> hazard:alpha1       1.032208   0.360466   2.864  0.00419 ** 
#> hazard:alpha2      -1.503094   0.317228  -4.738 2.16e-06 ***
#> hazard:phi1         0.774640   0.141422   5.478 4.31e-08 ***
#> hazard:phi2        -1.223616   0.143948  -8.500  < 2e-16 ***
#> longitudinal:beta1 -1.021807   0.014190 -72.009  < 2e-16 ***
#> longitudinal:beta2 -0.602169   0.019091 -31.543  < 2e-16 ***
#> longitudinal:beta3 -0.006041   0.006801  -0.888  0.37440    
#> longitudinal:beta4 -0.805953   0.012942 -62.273  < 2e-16 ***
#> longitudinal:beta5  0.507161   0.008583  59.087  < 2e-16 ***
#> longitudinal:beta6  0.406766   0.006644  61.225  < 2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> ---
#> Log-likelihood: 931.3818   AIC: -1832.764   BIC: -1793.686
#> N = 100  Convergence: EM algorithm converged after 6 iterations

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
