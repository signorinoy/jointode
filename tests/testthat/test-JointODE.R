test_that("JointODE basic functionality", {
  # Quick test that runs on CI/codecov for coverage
  # Using minimal data and iterations to keep test fast

  # Only skip on CRAN where time limits are strict
  skip_on_cran()

  # Generate minimal test data - very small for speed
  tiny_data <- JointODE:::.create_example_data(
    n_subjects = 50, # Minimum viable data
    seed = 123
  )

  # Test that function can be called and returns expected structure
  result <- JointODE(
    longitudinal_formula = observed ~ x1 + x2,
    survival_formula = Surv(time, status) ~ w1 + w2,
    longitudinal_data = tiny_data$data$longitudinal_data,
    survival_data = tiny_data$data$survival_data,
    state = as.matrix(tiny_data$data$state),
    control = list(em_maxit = 2)
  )

  # Basic structure checks
  expect_s3_class(result, "JointODE")
  expect_type(result, "list")

  # Check all required components exist based on actual return structure
  required_components <- c(
    "parameters",
    "logLik",
    "AIC",
    "BIC",
    "cindex",
    "convergence",
    "random_effects",
    "vcov",
    "data",
    "control",
    "robust",
    "call"
  )
  for (comp in required_components) {
    expect_true(
      comp %in% names(result),
      info = paste("Missing component:", comp)
    )
  }

  # Check parameters structure
  expect_type(result$parameters, "list")
  expect_true("coefficients" %in% names(result$parameters))
  expect_true("configurations" %in% names(result$parameters))

  # Check coefficients structure
  coef_names <- c(
    "baseline",
    "acceleration",
    "hazard",
    "measurement_error_sd",
    "random_effect_sd"
  )
  for (name in coef_names) {
    expect_true(
      name %in% names(result$parameters$coefficients),
      info = paste("Missing coefficient:", name)
    )
  }

  # Check vcov - may be NULL if not converged or robust estimation failed
  if (!is.null(result$vcov)) {
    n_params <- length(unlist(result$parameters$coefficients[c(
      "baseline",
      "acceleration",
      "hazard"
    )]))
    expect_true(is.matrix(result$vcov))
    expect_equal(dim(result$vcov), c(n_params, n_params))
  }

  # Check logLik, AIC, BIC
  expect_type(result$logLik, "double")
  expect_true(is.finite(result$logLik))
  expect_type(result$AIC, "double")
  expect_type(result$BIC, "double")

  # Check convergence info
  expect_type(result$convergence, "list")
  expect_true("converged" %in% names(result$convergence))
  expect_true("em_iterations" %in% names(result$convergence))
  expect_true("message" %in% names(result$convergence))

  # Check random effects structure
  expect_type(result$random_effects, "list")
  expect_true("estimates" %in% names(result$random_effects))
  expect_true("variances" %in% names(result$random_effects))

  # Check C-index if exists
  if (!is.null(result$cindex)) {
    expect_type(result$cindex, "double")
    expect_true(result$cindex >= 0 && result$cindex <= 1)
  }

  # Check methods work
  expect_output(print(result), "JointODE")
  expect_s3_class(summary(result), "summary.JointODE")
  expect_type(coef(result), "double")

  # vcov method may return NULL if vcov is NULL
  vcov_result <- vcov(result)
  if (!is.null(vcov_result)) {
    expect_true(is.matrix(vcov_result))
  }
})
