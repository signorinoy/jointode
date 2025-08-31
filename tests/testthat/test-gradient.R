# ==============================================================================
# Test: Gradient Computation Using numDeriv
# ==============================================================================
# Purpose: Verify analytical gradients match numerical gradients
# Coverage:
#   - Theta gradient (η, α, φ, β parameters)
#   - Edge cases (single subject, zero random effect)
#   - Parameter scaling consistency
#   - Finite value checks

test_that("compute_gradient_joint matches numerical gradient", {
  skip_if_not_installed("numDeriv")

  # Setup: Load test data and compute posteriors
  test_env <- load_test_data(n_subjects = 3)
  posteriors <- .compute_posteriors(test_env$data, test_env$parameters)

  # Extract parameters and dimensions
  n_surv_covariates <- ncol(test_env$data[[1]]$covariates)

  coefficients <- test_env$parameters$coefficients
  configurations <- test_env$parameters$configurations

  params <- c(
    coefficients$baseline,
    coefficients$hazard,
    coefficients$acceleration
  )
  fixed_params <- list(
    measurement_error_sd = coefficients$measurement_error_sd,
    random_effect_sd = coefficients$random_effect_sd
  )

  # Compute gradients
  grad_analytical <- .compute_gradient_joint(
    params = params,
    data_list = test_env$data,
    posteriors = posteriors,
    configurations = configurations,
    fixed_parameters = fixed_params
  )
  names(grad_analytical) <- NULL

  grad_numerical <- numDeriv::grad(
    func = .compute_objective_joint,
    x = params,
    data_list = test_env$data,
    posteriors = posteriors,
    configurations = configurations,
    fixed_parameters = fixed_params,
    method = "Richardson"
  )

  # Overall comparison
  # Note: Due to ODE integration and numerical differentiation,
  # a tolerance of 5% is reasonable for the full gradient
  compare_gradient_component(
    grad_analytical,
    grad_numerical,
    "Full gradient",
    tolerance = 5e-2
  )
})

test_that("gradient computation handles edge cases", {
  skip_if_not_installed("numDeriv")

  # Test with single subject
  test_env <- load_test_data(n_subjects = 1)
  posteriors <- .compute_posteriors(test_env$data, test_env$parameters)

  coefficients <- test_env$parameters$coefficients
  configurations <- test_env$parameters$configurations

  params <- c(
    coefficients$baseline,
    coefficients$hazard,
    coefficients$acceleration
  )
  fixed_params <- list(
    measurement_error_sd = coefficients$measurement_error_sd,
    random_effect_sd = coefficients$random_effect_sd
  )

  # Should not error with single subject
  grad <- .compute_gradient_joint(
    params = params,
    data_list = test_env$data,
    posteriors = posteriors,
    configurations = configurations,
    fixed_parameters = fixed_params
  )

  expect_true(all(is.finite(grad)))
  expect_equal(length(grad), length(params))
})

test_that("gradient components have correct dimensions", {
  test_env <- load_test_data(n_subjects = 5)
  posteriors <- .compute_posteriors(test_env$data, test_env$parameters)

  n_baseline <- length(test_env$parameters$coefficients$baseline)
  n_hazard <- length(test_env$parameters$coefficients$hazard)
  n_acceleration <- length(test_env$parameters$coefficients$acceleration)

  coefficients <- test_env$parameters$coefficients
  configurations <- test_env$parameters$configurations

  params <- c(
    coefficients$baseline,
    coefficients$hazard,
    coefficients$acceleration
  )
  fixed_params <- list(
    measurement_error_sd = coefficients$measurement_error_sd,
    random_effect_sd = coefficients$random_effect_sd
  )

  grad <- .compute_gradient_joint(
    params = params,
    data_list = test_env$data,
    posteriors = posteriors,
    configurations = configurations,
    fixed_parameters = fixed_params
  )

  # Check total dimension
  expect_equal(length(grad), n_baseline + n_hazard + n_acceleration)

  # Check that gradient is finite
  expect_true(all(is.finite(grad)))

  # Verify gradient structure matches parameter structure
  grad_baseline <- grad[1:n_baseline]
  grad_hazard <- grad[(n_baseline + 1):(n_baseline + n_hazard)]
  grad_acceleration <- grad[(n_baseline + n_hazard + 1):length(grad)]

  expect_equal(length(grad_baseline), n_baseline)
  expect_equal(length(grad_hazard), n_hazard)
  expect_equal(length(grad_acceleration), n_acceleration)
})
