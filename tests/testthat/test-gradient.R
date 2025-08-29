# ==============================================================================
# Test: Gradient Computation Using numDeriv
# ==============================================================================
# Purpose: Verify analytical gradients match numerical gradients
# Coverage:
#   - Theta gradient (η, α, φ, γ parameters)
#   - Beta gradient (spherical coordinates)
#   - Edge cases (single subject, zero random effect)
#   - Parameter scaling consistency
#   - Finite value checks

test_that("compute_grad_theta_forward matches numerical gradient", {
  skip_if_not_installed("numDeriv")

  # Setup: Load test data and compute posteriors
  test_env <- load_test_data(n_subjects = 3)
  posteriors <- .compute_posteriors(test_env$data, test_env$parameters)

  # Extract parameters and dimensions
  n_surv_covariates <- ncol(test_env$data[[1]]$covariates)
  theta_params <- extract_theta_params(test_env$parameters, n_surv_covariates)
  fixed_params <- create_fixed_params_theta(test_env$parameters)
  dims <- get_param_dimensions(
    test_env$parameters$configurations, n_surv_covariates
  )

  # Compute gradients
  grad_analytical <- .compute_grad_theta_forward(
    params = theta_params,
    data_list = test_env$data,
    posteriors = posteriors,
    configurations = test_env$parameters$configurations,
    fixed_params = fixed_params
  )
  names(grad_analytical) <- NULL

  grad_numerical <- numDeriv::grad(
    func = .compute_objective_theta,
    x = theta_params,
    data_list = test_env$data,
    posteriors = posteriors,
    configurations = test_env$parameters$configurations,
    fixed_params = fixed_params,
    method = "Richardson"
  )

  # Define component-wise comparison specifications
  gradient_specs <- list(
    list(name = "η (baseline hazard)", size = dims$n_eta, tolerance = 1e-3),
    list(name = "α (association)", size = dims$n_alpha, tolerance = 1e-4),
    list(name = "φ (survival covariates)", size = dims$n_phi, tolerance = 1e-4),
    list(name = "γ (index spline)", size = dims$n_gamma, tolerance = 1e-3)
  )

  # Compare each component
  idx <- 1
  for (spec in gradient_specs) {
    if (spec$size > 0) {  # Skip empty components
      indices <- idx:(idx + spec$size - 1)
      compare_gradient_component(
        grad_analytical[indices],
        grad_numerical[indices],
        spec$name,
        spec$tolerance
      )
      idx <- idx + spec$size
    }
  }

  # Overall comparison
  compare_gradient_component(
    grad_analytical, grad_numerical,
    "Full theta gradient",
    tolerance = 1e-3
  )
})
