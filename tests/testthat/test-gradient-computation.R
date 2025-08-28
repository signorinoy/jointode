# Tests for gradient computation using numDeriv
# Verify analytical gradients against numerical gradients

test_that("compute_grad_theta_forward matches numerical gradient", {
  skip_if_not_installed("numDeriv")

  # Setup test environment using helper functions
  test_env <- setup_test_environment()
  data_list <- test_env$data_list
  config <- test_env$config
  posteriors <- test_env$posteriors

  # Take only first few subjects for faster testing
  n_test_subjects <- min(3, length(data_list))
  data_list_subset <- data_list[1:n_test_subjects]
  posteriors_subset <- list(
    b = posteriors$b[1:n_test_subjects],
    v = posteriors$v[1:n_test_subjects],
    exp_b = posteriors$exp_b[1:n_test_subjects]
  )

  # Fixed parameters for theta optimization
  fixed_params_theta <- list(
    index_beta = test_env$params$index_beta,
    measurement_error_sd = test_env$params$measurement_error_sd,
    random_effect_sd = test_env$params$random_effect_sd
  )

  # Get number of survival covariates
  n_surv_covariates <- ncol(data_list_subset[[1]]$covariates)

  # Create theta parameter vector (η, α, φ, γ)
  n_eta <- config$baseline$df
  n_alpha <- 3
  n_phi <- n_surv_covariates
  n_gamma <- config$index$df

  theta_params <- c(
    test_env$params$baseline,      # η
    test_env$params$hazard[1:3],   # α
    test_env$params$hazard[4:(3 + n_phi)], # φ
    test_env$params$index_g        # γ
  )

  # Compute analytical gradient
  grad_analytical <- .compute_grad_theta_forward(
    params = theta_params,
    data_list = data_list_subset,
    posteriors = posteriors_subset,
    config = config,
    fixed_params = fixed_params_theta
  )
  names(grad_analytical) <- NULL

  # Compute numerical gradient using numDeriv
  grad_numerical <- numDeriv::grad(
    func = .compute_objective_theta,
    x = theta_params,
    data_list = data_list_subset,
    posteriors = posteriors_subset,
    config = config,
    fixed_params = fixed_params_theta,
    method = "Richardson"
  )

  # Compare gradients component by component for better diagnostics
  idx <- 1

  # Test η (baseline hazard parameters)
  eta_analytical <- grad_analytical[idx:(idx + n_eta - 1)]
  eta_numerical <- grad_numerical[idx:(idx + n_eta - 1)]
  expect_equal(eta_analytical, eta_numerical,
               tolerance = 1e-3,
               label = "η gradient",
               info = paste("Max difference:",
                            max(abs(eta_analytical - eta_numerical))))
  idx <- idx + n_eta

  # Test α (association parameters)
  alpha_analytical <- grad_analytical[idx:(idx + n_alpha - 1)]
  alpha_numerical <- grad_numerical[idx:(idx + n_alpha - 1)]
  expect_equal(alpha_analytical, alpha_numerical,
               tolerance = 1e-4,
               label = "α gradient",
               info = paste("Max difference:",
                            max(abs(alpha_analytical - alpha_numerical))))
  idx <- idx + n_alpha

  # Test φ (survival covariates) if present
  if (n_phi > 0) {
    phi_analytical <- grad_analytical[idx:(idx + n_phi - 1)]
    phi_numerical <- grad_numerical[idx:(idx + n_phi - 1)]
    expect_equal(phi_analytical, phi_numerical,
                 tolerance = 1e-4,
                 label = "φ gradient",
                 info = paste("Max difference:",
                              max(abs(phi_analytical - phi_numerical))))
    idx <- idx + n_phi
  }

  # Test γ (index spline parameters)
  gamma_analytical <- grad_analytical[idx:(idx + n_gamma - 1)]
  gamma_numerical <- grad_numerical[idx:(idx + n_gamma - 1)]
  expect_equal(gamma_analytical, gamma_numerical,
               tolerance = 1e-3,  # Tighter tolerance with proper scaling
               label = "γ gradient",
               info = paste("Max difference:",
                            max(abs(gamma_analytical - gamma_numerical))))

  # Overall gradient comparison
  expect_equal(grad_analytical, grad_numerical,
               tolerance = 1e-3,  # Tighter tolerance with proper scaling
               label = "Full theta gradient",
               info = paste("Max difference:",
                            max(abs(grad_analytical - grad_numerical)),
                            "\nRelative error:",
                            max(abs((grad_analytical - grad_numerical) /
                                      (abs(grad_numerical) + 1e-10)))))
})

test_that("compute_grad_beta_forward matches numerical gradient", {
  skip_if_not_installed("numDeriv")

  # Setup test environment
  test_env <- setup_test_environment()
  data_list <- test_env$data_list
  config <- test_env$config
  posteriors <- test_env$posteriors

  # Take only first few subjects for faster testing
  n_test_subjects <- min(3, length(data_list))
  data_list_subset <- data_list[1:n_test_subjects]
  posteriors_subset <- list(
    b = posteriors$b[1:n_test_subjects],
    v = posteriors$v[1:n_test_subjects],
    exp_b = posteriors$exp_b[1:n_test_subjects]
  )

  # Fixed parameters for beta optimization
  fixed_params_beta <- list(
    baseline = test_env$params$baseline,
    hazard = test_env$params$hazard,
    index_g = test_env$params$index_g,
    measurement_error_sd = test_env$params$measurement_error_sd,
    random_effect_sd = test_env$params$random_effect_sd
  )

  # Convert beta to spherical coordinates
  beta_original <- test_env$params$index_beta
  spherical_params <- .beta_to_spherical(beta_original)

  # Compute analytical gradient
  grad_analytical <- .compute_grad_beta_forward(
    params = spherical_params,
    data_list = data_list_subset,
    posteriors = posteriors_subset,
    config = config,
    fixed_params = fixed_params_beta
  )

  # Compute numerical gradient using numDeriv
  grad_numerical <- numDeriv::grad(
    func = .compute_objective_beta,
    x = spherical_params,
    data_list = data_list_subset,
    posteriors = posteriors_subset,
    config = config,
    fixed_params = fixed_params_beta,
    method = "Richardson"
  )

  # Compare gradients
  expect_equal(grad_analytical, grad_numerical,
               tolerance = 1e-4,
               label = "Beta gradient (spherical coordinates)",
               info = paste("Max difference:",
                            max(abs(grad_analytical - grad_numerical)),
                            "\nRelative error:",
                            max(abs((grad_analytical - grad_numerical) /
                                      (abs(grad_numerical) + 1e-10)))))

  # Test that spherical to beta conversion and back works correctly
  beta_reconstructed <- .spherical_to_beta(spherical_params)
  expect_equal(beta_reconstructed, beta_original,
               tolerance = 1e-10,
               label = "Beta reconstruction from spherical")

  # Verify beta is unit vector
  expect_equal(sum(beta_reconstructed^2), 1,
               tolerance = 1e-10,
               label = "Beta unit vector constraint")
})

test_that("gradient computation handles edge cases correctly", {
  skip_if_not_installed("numDeriv")

  # Setup minimal test environment
  test_env <- setup_test_environment()
  data_list <- test_env$data_list[1] # Single subject
  config <- test_env$config
  posteriors <- list(
    b = test_env$posteriors$b[1],
    v = test_env$posteriors$v[1],
    exp_b = test_env$posteriors$exp_b[1]
  )

  # Test with single subject
  fixed_params_theta <- list(
    index_beta = test_env$params$index_beta,
    measurement_error_sd = test_env$params$measurement_error_sd,
    random_effect_sd = test_env$params$random_effect_sd
  )

  n_surv_covariates <- ncol(data_list[[1]]$covariates)
  n_eta <- config$baseline$df
  n_alpha <- 3
  n_phi <- n_surv_covariates
  n_gamma <- config$index$df

  theta_params <- c(
    test_env$params$baseline,
    test_env$params$hazard[1:3],
    test_env$params$hazard[4:(3 + n_phi)],
    test_env$params$index_g
  )

  # Compute gradients for single subject
  grad_analytical <- .compute_grad_theta_forward(
    params = theta_params,
    data_list = data_list,
    posteriors = posteriors,
    config = config,
    fixed_params = fixed_params_theta
  )

  grad_numerical <- numDeriv::grad(
    func = .compute_objective_theta,
    x = theta_params,
    data_list = data_list,
    posteriors = posteriors,
    config = config,
    fixed_params = fixed_params_theta,
    method = "Richardson"
  )

  expect_equal(grad_analytical, grad_numerical,
               tolerance = 1e-3,  # Slightly relaxed for edge cases
               label = "Single subject gradient")

  # Test with zero random effect
  posteriors_zero_b <- list(
    b = 0,
    v = posteriors$v,
    exp_b = 1
  )

  grad_analytical_zero_b <- .compute_grad_theta_forward(
    params = theta_params,
    data_list = data_list,
    posteriors = posteriors_zero_b,
    config = config,
    fixed_params = fixed_params_theta
  )

  grad_numerical_zero_b <- numDeriv::grad(
    func = .compute_objective_theta,
    x = theta_params,
    data_list = data_list,
    posteriors = posteriors_zero_b,
    config = config,
    fixed_params = fixed_params_theta,
    method = "Richardson"
  )

  expect_equal(grad_analytical_zero_b, grad_numerical_zero_b,
               tolerance = 1e-3,  # Slightly relaxed for edge cases
               label = "Zero random effect gradient")
})

test_that("gradient computation consistent across parameter scales", {
  skip_if_not_installed("numDeriv")

  # Setup test environment
  test_env <- setup_test_environment()
  data_list <- test_env$data_list[1:2] # Two subjects for speed
  config <- test_env$config
  posteriors <- list(
    b = test_env$posteriors$b[1:2],
    v = test_env$posteriors$v[1:2],
    exp_b = test_env$posteriors$exp_b[1:2]
  )

  fixed_params_theta <- list(
    index_beta = test_env$params$index_beta,
    measurement_error_sd = test_env$params$measurement_error_sd,
    random_effect_sd = test_env$params$random_effect_sd
  )

  n_surv_covariates <- ncol(data_list[[1]]$covariates)
  n_eta <- config$baseline$df
  n_alpha <- 3
  n_phi <- n_surv_covariates
  n_gamma <- config$index$df

  # Test with scaled parameters
  scales <- c(0.1, 1.0, 10.0)

  for (scale in scales) {
    theta_params <- c(
      test_env$params$baseline * scale,
      test_env$params$hazard[1:3] * scale,
      test_env$params$hazard[4:(3 + n_phi)] * scale,
      test_env$params$index_g * scale
    )

    grad_analytical <- .compute_grad_theta_forward(
      params = theta_params,
      data_list = data_list,
      posteriors = posteriors,
      config = config,
      fixed_params = fixed_params_theta
    )

    grad_numerical <- numDeriv::grad(
      func = .compute_objective_theta,
      x = theta_params,
      data_list = data_list,
      posteriors = posteriors,
      config = config,
      fixed_params = fixed_params_theta,
      method = "Richardson"
    )

    expect_equal(grad_analytical, grad_numerical,
                 tolerance = 5e-3,  # Relaxed for parameter scaling tests
                 label = paste("Gradient with scale", scale),
                 info = paste("Scale:", scale,
                              "Max difference:",
                              max(abs(grad_analytical - grad_numerical))))
  }
})

test_that("objective function and gradient are finite for reasonable inputs", {
  # Setup test environment
  test_env <- setup_test_environment()
  data_list <- test_env$data_list[1:2]
  config <- test_env$config
  posteriors <- list(
    b = test_env$posteriors$b[1:2],
    v = test_env$posteriors$v[1:2],
    exp_b = test_env$posteriors$exp_b[1:2]
  )

  fixed_params_theta <- list(
    index_beta = test_env$params$index_beta,
    measurement_error_sd = test_env$params$measurement_error_sd,
    random_effect_sd = test_env$params$random_effect_sd
  )

  n_surv_covariates <- ncol(data_list[[1]]$covariates)
  n_eta <- config$baseline$df
  n_alpha <- 3
  n_phi <- n_surv_covariates
  n_gamma <- config$index$df

  theta_params <- c(
    test_env$params$baseline,
    test_env$params$hazard[1:3],
    test_env$params$hazard[4:(3 + n_phi)],
    test_env$params$index_g
  )

  # Check objective function is finite
  obj_value <- .compute_objective_theta(
    params = theta_params,
    data_list = data_list,
    posteriors = posteriors,
    config = config,
    fixed_params = fixed_params_theta
  )

  expect_true(is.finite(obj_value),
              label = "Objective function is finite")

  # Check gradient is finite
  grad_value <- .compute_grad_theta_forward(
    params = theta_params,
    data_list = data_list,
    posteriors = posteriors,
    config = config,
    fixed_params = fixed_params_theta
  )

  expect_true(all(is.finite(grad_value)),
              label = "Gradient components are finite")

  # Test beta objective and gradient
  fixed_params_beta <- list(
    baseline = test_env$params$baseline,
    hazard = test_env$params$hazard,
    index_g = test_env$params$index_g,
    measurement_error_sd = test_env$params$measurement_error_sd,
    random_effect_sd = test_env$params$random_effect_sd
  )

  spherical_params <- .beta_to_spherical(test_env$params$index_beta)

  obj_beta <- .compute_objective_beta(
    params = spherical_params,
    data_list = data_list,
    posteriors = posteriors,
    config = config,
    fixed_params = fixed_params_beta
  )

  expect_true(is.finite(obj_beta),
              label = "Beta objective function is finite")

  grad_beta <- .compute_grad_beta_forward(
    params = spherical_params,
    data_list = data_list,
    posteriors = posteriors,
    config = config,
    fixed_params = fixed_params_beta
  )

  expect_true(all(is.finite(grad_beta)),
              label = "Beta gradient components are finite")
})
