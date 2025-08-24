# Tests for .solve_joint_ode function

test_that(".solve_joint_ode works with basic sensitivity type", {
  # Create test data structure
  test_data <- simulate()
  test_data_process <- .process(
    longitudinal_data = test_data$longitudinal_data,
    longitudinal_formula = v ~ x1 + x2,
    survival_data = test_data$survival_data,
    survival_formula = Surv(time, status) ~ w1 + w2,
    id = "id",
    time = "time"
  )
  n_subjects <- length(test_data_process)

  # Create test parameters
  g_0 <- function(u) 2 * tanh(u / 3)
  lambda_0 <- function(t) log((1.5 / 8) * (t / 8)^0.5)

  times <- test_data$survival_data[, "time"]
  spline_baseline <- list(
    degree = 3,
    n_knots = 5,
    knot_placement = "quantile",
    boundary_knots = NULL
  )
  baseline_spline_coefficients <- estimate_bspline_coef(
    times, lambda_0, spline_baseline
  )
  hazard_coefficients <- c(0.5, 0.4, -0.5, 0.2, -0.15)

  index_coefficients <- c(-0.3, -0.5, 0, 0.2, 0.1, 0.05)
  Z <- test_data$longitudinal_data[
    , c("biomarker", "velocity", "x1", "x2", "time")
  ]
  u <- as.matrix(Z) %*% as.matrix(c(-0.3, -0.5, 0.2, 0.1, 0.05))
  spline_index <- list(
    degree = 3,
    n_knots = 10,
    knot_placement = "quantile",
    boundary_knots = NULL
  )
  index_spline_coefficients <- estimate_bspline_coef(
    u, g_0, spline_index
  )
  index_coefficients <- index_coefficients / sqrt(sum(index_coefficients^2))
  measurement_error_sd <- 0.1
  random_effect_sd <- 0.1

  spline_baseline_config <- .get_spline_config(
    x = times,
    degree = spline_baseline$degree,
    n_knots = spline_baseline$n_knots,
    knot_placement = spline_baseline$knot_placement,
    boundary_knots = spline_baseline$boundary_knots
  )

  spline_index_config <- .get_spline_config(
    x = u,
    degree = spline_index$degree,
    n_knots = spline_index$n_knots,
    knot_placement = spline_index$knot_placement,
    boundary_knots = spline_index$boundary_knots
  )

  test_params <- list(
    coef = list(
      baseline = baseline_spline_coefficients,
      hazard = hazard_coefficients,
      index_g = index_spline_coefficients,
      index_beta = index_coefficients
    ),
    config = list(
      baseline = spline_baseline_config,
      index = spline_index_config
    )
  )

  # Test basic sensitivity type (default)
  for (i in seq_len(n_subjects)) {
    subject_data <- test_data_process[[i]]
    subject_id <- names(test_data_process)[i]
    ode_solution <- .solve_joint_ode(subject_data, test_params)
    idx <- test_data$longitudinal_data$id == subject_id
    biomarker_data <- test_data$longitudinal_data[idx, ]
    biomarker_true <- biomarker_data$biomarker
    velocity_true <- biomarker_data$velocity
    acceleration_true <- biomarker_data$acceleration
    expect_equal(ode_solution$biomarker, biomarker_true, tolerance = 1e-2)
    expect_equal(ode_solution$velocity, velocity_true, tolerance = 1e-3)
    expect_equal(ode_solution$acceleration, acceleration_true, tolerance = 1e-3)
  }
})
