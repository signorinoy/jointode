# ==============================================================================
# Tests for Simulation Functions
# (simulate, .estimate_bspline_coef, .create_example_data)
# ==============================================================================
#
# Test Coverage:
#
# 1. simulate() Function Tests
#    - Data structure validation (longitudinal and survival components)
#    - Parameter handling and validation
#    - Input validation (negative values, wrong dimensions)
#    - Edge cases (small/large sample sizes)
#    - Reasonable value ranges
#    - Reproducibility with seeds
#    - Verbose option functionality
#
# 2. .estimate_bspline_coef() Function Tests
#    - Linear function approximation accuracy
#    - Nonlinear function approximation (sine, polynomial)
#    - Constant function handling
#    - Small sample size robustness
#    - Coefficient properties validation
#
# 3. .create_example_data() Function Tests
#    - Complete data structure generation
#    - Formula inclusion
#    - Parameter structure validation
#    - Configuration correctness
#    - Different sample size handling
#    - Reproducibility with seeds
#    - Parameter normalization verification
# ==============================================================================

# Tests for main simulate() function
test_that("simulate() generates valid data structure", {
  set.seed(42)
  sim_data <- simulate(n = 30)

  # Check overall structure
  expect_type(sim_data, "list")
  expect_named(sim_data, c("longitudinal_data", "survival_data"))

  # Check longitudinal data
  long_data <- sim_data$longitudinal_data
  expect_s3_class(long_data, "data.frame")
  expect_true(all(
    c(
      "id",
      "time",
      "v",
      "x1",
      "x2",
      "biomarker",
      "velocity",
      "acceleration"
    ) %in%
      names(long_data)
  ))

  # Check data types
  expect_type(long_data$id, "integer")
  expect_type(long_data$time, "double")
  expect_type(long_data$v, "double")
  expect_type(long_data$x1, "double")
  expect_type(long_data$x2, "double")

  # Check survival data
  surv_data <- sim_data$survival_data
  expect_s3_class(surv_data, "data.frame")
  expect_true(all(
    c("id", "time", "status", "w1", "w2", "b") %in%
      names(surv_data)
  ))
  expect_equal(nrow(surv_data), 30)

  # Check status values
  expect_true(all(surv_data$status %in% c(0, 1)))
})

test_that("simulate() respects input parameters", {
  set.seed(123)

  # Test with custom parameters
  custom_alpha <- c(0.3, 0.2)
  custom_beta <- c(-0.2, -0.3, 0.1, 0.15, 0.08)
  custom_phi <- c(0.3, -0.2)

  sim_data <- simulate(
    n = 25,
    alpha = custom_alpha,
    beta = custom_beta,
    phi = custom_phi,
    weibull_shape = 2,
    weibull_scale = 10,
    sigma_b = 0.15,
    sigma_e = 0.08,
    seed = 456,
    verbose = FALSE
  )

  # Check sample size
  expect_equal(nrow(sim_data$survival_data), 25)
  expect_equal(length(unique(sim_data$longitudinal_data$id)), 25)

  # Check that parameters affect the data
  # (different seed should give different results)
  sim_data2 <- simulate(n = 25, seed = 789)
  expect_false(all(sim_data$survival_data$time == sim_data2$survival_data$time))
})

test_that("simulate() handles edge cases", {
  # Very small sample size
  set.seed(111)
  small_sim <- simulate(n = 2)
  expect_equal(nrow(small_sim$survival_data), 2)

  # Small sample size
  set.seed(112)
  small_sim2 <- simulate(n = 5)
  expect_equal(nrow(small_sim2$survival_data), 5)

  # Moderate sample size
  set.seed(222)
  large_sim <- simulate(n = 50)
  expect_equal(nrow(large_sim$survival_data), 50)

  # Check all subjects have longitudinal data
  surv_ids <- large_sim$survival_data$id
  long_ids <- unique(large_sim$longitudinal_data$id)
  expect_true(all(surv_ids %in% long_ids))
})

test_that("simulate() produces reasonable values", {
  set.seed(333)
  sim_data <- simulate(n = 50)

  # Check for reasonable survival times
  expect_true(all(sim_data$survival_data$time > 0))
  expect_true(all(sim_data$survival_data$time < 100)) # Should not be extreme

  # Check for reasonable event rate
  # (should be between 20% and 90% for typical settings)
  event_rate <- mean(sim_data$survival_data$status)
  expect_true(event_rate > 0.2 && event_rate < 0.9)

  # Check longitudinal measurements are reasonable
  expect_false(any(is.na(sim_data$longitudinal_data$v)))
  expect_false(any(is.infinite(sim_data$longitudinal_data$v)))

  # Check that velocity and acceleration are computed
  expect_false(all(sim_data$longitudinal_data$velocity == 0))
  expect_false(all(sim_data$longitudinal_data$acceleration == 0))
})

test_that("simulate() input validation works", {
  # Invalid n
  expect_error(simulate(n = -5), "n must be a positive integer")
  expect_error(simulate(n = 0), "n must be a positive integer")
  expect_error(simulate(n = 2.5), "n must be a positive integer")

  # Invalid alpha
  expect_error(
    simulate(n = 10, alpha = c(0.5)),
    "alpha must be a numeric vector of length 2"
  )

  # Invalid beta
  expect_error(
    simulate(n = 10, beta = c(-0.3, -0.5)),
    "beta must be a numeric vector of length 5"
  )

  # Invalid phi
  expect_error(
    simulate(n = 10, phi = c(0.2)),
    "phi must be a numeric vector of length 2"
  )

  # Invalid weibull parameters
  expect_error(
    simulate(n = 10, weibull_shape = -1),
    "weibull_shape must be positive"
  )
  expect_error(
    simulate(n = 10, weibull_scale = 0),
    "weibull_scale must be positive"
  )

  # Invalid standard deviations
  expect_error(simulate(n = 10, sigma_b = -0.1), "sigma_b must be positive")
  expect_error(simulate(n = 10, sigma_e = 0), "sigma_e must be positive")
})

test_that("simulate() verbose option works", {
  # Test verbose output
  expect_message(
    simulate(n = 10, verbose = TRUE),
    "Step 1/2: Generating survival times"
  )

  # Should also show Step 2
  expect_message(
    simulate(n = 10, verbose = TRUE),
    "Step 2/2: Generating longitudinal data"
  )

  # Test silent mode
  expect_silent(simulate(n = 10))
})

test_that("simulate() produces consistent results with seed", {
  # Same seed should produce identical results
  set.seed(999)
  sim1 <- simulate(n = 20)

  set.seed(999)
  sim2 <- simulate(n = 20)

  # Check longitudinal data
  expect_equal(sim1$longitudinal_data$v, sim2$longitudinal_data$v)

  # Check survival data
  expect_equal(sim1$survival_data$time, sim2$survival_data$time)
  expect_equal(sim1$survival_data$status, sim2$survival_data$status)
  expect_equal(sim1$survival_data$b, sim2$survival_data$b)
})

# Tests for helper functions

test_that(".estimate_bspline_coef approximates functions correctly", {
  # Test 1: Linear function approximation
  x <- seq(0, 10, length.out = 50)
  f_linear <- function(t) 2 * t + 3
  config <- list(
    degree = 3,
    n_knots = 5,
    knot_placement = "quantile",
    boundary_knots = NULL
  )

  coef <- .estimate_bspline_coef(x, f_linear, config)

  # Coefficients should be numeric
  expect_type(coef, "double")
  expect_true(length(coef) > 0)

  # Reconstruct and check approximation quality
  spline_config <- .get_spline_config(
    x = x,
    degree = config$degree,
    n_knots = config$n_knots,
    knot_placement = config$knot_placement,
    boundary_knots = config$boundary_knots
  )
  basis <- .compute_spline_basis(x, spline_config)
  y_approx <- as.vector(basis %*% coef)
  y_true <- f_linear(x)

  # Should approximate well (R-squared > 0.99 for linear function)
  r_squared <- 1 - sum((y_true - y_approx)^2) / sum((y_true - mean(y_true))^2)
  expect_gt(r_squared, 0.99)
})

test_that(".estimate_bspline_coef handles nonlinear functions", {
  # Test 2: Sine function approximation
  x <- seq(0, 2 * pi, length.out = 100)
  f_sine <- function(t) sin(t)
  config <- list(
    degree = 3,
    n_knots = 10,
    knot_placement = "quantile",
    boundary_knots = NULL
  )

  coef <- .estimate_bspline_coef(x, f_sine, config)

  # Check coefficient properties
  expect_type(coef, "double")
  expect_equal(length(coef), config$n_knots + config$degree + 1)
  expect_false(any(is.na(coef)))
  expect_false(any(is.infinite(coef)))

  # Verify approximation quality
  spline_config <- .get_spline_config(
    x = x,
    degree = config$degree,
    n_knots = config$n_knots,
    knot_placement = config$knot_placement,
    boundary_knots = config$boundary_knots
  )
  basis <- .compute_spline_basis(x, spline_config)
  y_approx <- as.vector(basis %*% coef)
  y_true <- f_sine(x)

  # Should have reasonable approximation (R-squared > 0.95)
  r_squared <- 1 - sum((y_true - y_approx)^2) / sum((y_true - mean(y_true))^2)
  expect_gt(r_squared, 0.95)
})

test_that(".estimate_bspline_coef handles edge cases", {
  # Test 3: Constant function
  x <- seq(1, 5, length.out = 20)
  f_const <- function(t) rep(5, length(t)) # Ensure it returns a vector
  config <- list(
    degree = 2,
    n_knots = 3,
    knot_placement = "equal",
    boundary_knots = NULL
  )

  coef <- .estimate_bspline_coef(x, f_const, config)

  # Should handle constant function
  expect_type(coef, "double")
  expect_false(any(is.na(coef)))

  # Test 4: Small sample size
  x_small <- c(1, 2, 3, 4, 5)
  f_simple <- function(t) t^2
  config_small <- list(
    degree = 2,
    n_knots = 1,
    knot_placement = "quantile",
    boundary_knots = NULL
  )

  coef_small <- .estimate_bspline_coef(x_small, f_simple, config_small)
  expect_type(coef_small, "double")
  expect_true(length(coef_small) > 0)
})

test_that(".create_example_data generates valid dataset structure", {
  # Test with default parameters
  set.seed(123)
  example_data <- .create_example_data(n = 30)

  # Check overall structure
  expect_type(example_data, "list")
  expect_named(example_data, c("data", "formulas", "parameters"))

  # Check data component
  expect_type(example_data$data, "list")
  expect_named(example_data$data, c("longitudinal_data", "survival_data"))

  # Check longitudinal data
  long_data <- example_data$data$longitudinal_data
  expect_s3_class(long_data, "data.frame")
  expect_true(all(
    c(
      "id",
      "time",
      "v",
      "x1",
      "x2",
      "biomarker",
      "velocity",
      "acceleration"
    ) %in%
      names(long_data)
  ))

  # Check survival data
  surv_data <- example_data$data$survival_data
  expect_s3_class(surv_data, "data.frame")
  expect_true(all(
    c("id", "time", "status", "w1", "w2", "b") %in%
      names(surv_data)
  ))
  expect_equal(nrow(surv_data), 30)

  # Check formulas component
  expect_type(example_data$formulas, "list")
  expect_named(example_data$formulas, c("longitudinal", "survival"))
  expect_s3_class(example_data$formulas$longitudinal, "formula")
  expect_s3_class(example_data$formulas$survival, "formula")

  # Check parameters component
  expect_type(example_data$parameters, "list")
  expect_named(example_data$parameters, c("coefficients", "configurations"))

  # Check coefficients
  coef <- example_data$parameters$coefficients
  expect_type(coef, "list")
  expect_true(all(
    c(
      "baseline",
      "hazard",
      "acceleration",
      "measurement_error_sd",
      "random_effect_sd"
    ) %in%
      names(coef)
  ))

  # Check specific coefficient properties
  expect_type(coef$hazard, "double")
  expect_equal(length(coef$hazard), 4)
  expect_equal(coef$measurement_error_sd, 0.1)
  expect_equal(coef$random_effect_sd, 0.1)

  # Check configurations
  config <- example_data$parameters$configurations
  expect_type(config, "list")
  expect_named(config, "baseline")
  expect_type(config$baseline, "list")
})

test_that(".create_example_data handles different sample sizes", {
  # Small sample
  set.seed(456)
  small_data <- .create_example_data(n = 10)
  expect_equal(nrow(small_data$data$survival_data), 10)
  expect_equal(length(unique(small_data$data$longitudinal_data$id)), 10)

  # Larger sample
  set.seed(789)
  large_data <- .create_example_data(n = 100)
  expect_equal(nrow(large_data$data$survival_data), 100)
  expect_equal(length(unique(large_data$data$longitudinal_data$id)), 100)

  # Check that all subjects have longitudinal measurements
  long_ids <- unique(large_data$data$longitudinal_data$id)
  surv_ids <- large_data$data$survival_data$id
  expect_true(all(surv_ids %in% long_ids))
})

test_that(".create_example_data produces consistent results with seed", {
  # Generate data twice with same seed
  set.seed(999)
  data1 <- .create_example_data(n = 20)

  set.seed(999)
  data2 <- .create_example_data(n = 20)

  # Should be identical
  expect_equal(data1$data$survival_data$time, data2$data$survival_data$time)
  expect_equal(data1$data$survival_data$status, data2$data$survival_data$status)
  expect_equal(
    data1$parameters$coefficients$hazard,
    data2$parameters$coefficients$hazard
  )
})

test_that(".create_example_data validates parameter normalization", {
  example_data <- .create_example_data(n = 25)

  # Check acceleration coefficients
  acceleration <- example_data$parameters$coefficients$acceleration
  expect_true(all(is.finite(acceleration)))

  # Check B-spline coefficients exist and are finite
  baseline_coef <- example_data$parameters$coefficients$baseline
  expect_true(all(is.finite(baseline_coef)))
  expect_true(length(baseline_coef) > 0)
})
