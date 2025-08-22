# Tests for .solve_joint_ode function

test_that(".solve_joint_ode works with basic sensitivity type", {
  # Create test data structure
  test_data <- list(
    time = 5,
    status = 1,
    longitudinal = list(
      times = c(0, 1, 2, 3),
      measurements = c(1.0, 1.2, 1.5, 1.8),
      n_obs = 4,
      covariates = matrix(c(0.5, 0.6, 0.7, 0.8), ncol = 1)
    ),
    covariates = matrix(c(0.5, 1.2), ncol = 2)
  )

  # Create test parameters
  test_params <- list(
    coef = list(
      baseline = rep(0.1, 7),
      hazard = c(0.2, 0.1, 0.05, 0.3, 0.4),
      index_g = rep(0.05, 7),
      index_beta = c(0.01, 0.02, 0.03, 0.04)  # Small values
    ),
    config = list(
      baseline = list(
        degree = 3,
        knots = c(1, 2, 3),
        boundary_knots = c(0, 6),
        df = 7
      ),
      index = list(
        degree = 3,
        knots = c(-1, 0, 1),
        boundary_knots = c(-5, 5),
        df = 7
      )
    )
  )

  # Test basic sensitivity type (default)
  result <- JointODE:::.solve_joint_ode(test_data, test_params)

  expect_type(result, "list")
  expect_named(result, c("cum_hazard", "log_hazard", "biomarker"))
  expect_type(result$cum_hazard, "double")
  expect_type(result$log_hazard, "double")
  expect_length(result$biomarker, test_data$longitudinal$n_obs)
  expect_true(result$cum_hazard >= 0)

  # Test explicit basic sensitivity type
  result2 <- JointODE:::.solve_joint_ode(test_data, test_params,
                                         sensitivity_type = "basic")

  expect_equal(result$cum_hazard, result2$cum_hazard)
  expect_equal(result$log_hazard, result2$log_hazard)
  expect_equal(result$biomarker, result2$biomarker)
})

test_that(".solve_joint_ode validates sensitivity_type", {
  test_data <- list(
    time = 3,
    status = 1,
    longitudinal = list(
      times = c(0, 1),
      measurements = c(1.0, 1.2),
      n_obs = 2,
      covariates = NULL
    ),
    covariates = NULL
  )

  test_params <- list(
    coef = list(
      baseline = rep(0.1, 4),
      hazard = c(0.1, 0.05, 0.02),
      index_g = rep(0.05, 4),
      index_beta = c(0.01, 0.02, 0.005)  # Small values
    ),
    config = list(
      baseline = list(
        degree = 2,
        knots = c(1),
        boundary_knots = c(0, 4),
        df = 4
      ),
      index = list(
        degree = 2,
        knots = c(0),
        boundary_knots = c(-1, 1),
        df = 4
      )
    )
  )

  # Test invalid sensitivity type
  expect_error(
    JointODE:::.solve_joint_ode(test_data, test_params,
                                sensitivity_type = "invalid"),
    "Invalid sensitivity_type"
  )

  # Test valid types don't error
  expect_no_error(
    JointODE:::.solve_joint_ode(test_data, test_params,
                                sensitivity_type = "basic")
  )
})

test_that(".solve_joint_ode produces consistent results", {
  # Create reproducible test case
  set.seed(123)
  test_data <- list(
    time = 3.5,
    status = 1,
    longitudinal = list(
      times = seq(0, 3, by = 0.5),
      measurements = c(1.0, 1.1, 1.3, 1.4, 1.6, 1.8, 2.0),
      n_obs = 7,
      covariates = matrix(rnorm(7), ncol = 1)
    ),
    covariates = matrix(c(0.5), ncol = 1)
  )

  test_params <- list(
    coef = list(
      baseline = c(0.1, 0.12, 0.08, 0.15, 0.11, 0.09),
      hazard = c(0.2, 0.1, 0.05, 0.3),
      index_g = c(0.05, 0.06, 0.04, 0.07, 0.03, 0.02),
      index_beta = c(0.01, 0.02, 0.015, 0.025)  # Small values
    ),
    config = list(
      baseline = list(
        degree = 3,
        knots = c(1, 2),
        boundary_knots = c(0, 4),
        df = 6
      ),
      index = list(
        degree = 3,
        knots = c(-0.5, 0.5),
        boundary_knots = c(-5, 5),
        df = 6
      )
    )
  )

  # Run twice and check consistency
  result1 <- JointODE:::.solve_joint_ode(test_data, test_params, "basic")
  result2 <- JointODE:::.solve_joint_ode(test_data, test_params, "basic")

  expect_equal(result1$cum_hazard, result2$cum_hazard, tolerance = 1e-10)
  expect_equal(result1$log_hazard, result2$log_hazard, tolerance = 1e-10)
  expect_equal(result1$biomarker, result2$biomarker, tolerance = 1e-10)
})

test_that(".solve_joint_ode handles different event times correctly", {
  # Test with event time before last measurement
  early_event_data <- list(
    time = 1.5,
    status = 1,
    longitudinal = list(
      times = c(0, 1, 2, 3),
      measurements = c(1.0, 1.2, 1.4, 1.6),
      n_obs = 4,
      covariates = NULL
    ),
    covariates = NULL
  )

  # Test with event time after last measurement
  late_event_data <- list(
    time = 5,
    status = 0,
    longitudinal = list(
      times = c(0, 1, 2),
      measurements = c(1.0, 1.2, 1.4),
      n_obs = 3,
      covariates = NULL
    ),
    covariates = NULL
  )

  test_params <- list(
    coef = list(
      baseline = rep(0.1, 7),
      hazard = c(0.1, 0.05, 0.02),
      index_g = rep(0.05, 7),
      index_beta = c(0.01, 0.02, 0.005)  # Small values
    ),
    config = list(
      baseline = list(
        degree = 3,
        knots = c(1, 2, 3),
        boundary_knots = c(0, 6),
        df = 7
      ),
      index = list(
        degree = 3,
        knots = c(-1, 0, 1),
        boundary_knots = c(-5, 5),
        df = 7
      )
    )
  )

  # Both should work without errors
  result_early <- JointODE:::.solve_joint_ode(early_event_data, test_params)
  result_late <- JointODE:::.solve_joint_ode(late_event_data, test_params)

  expect_type(result_early, "list")
  expect_type(result_late, "list")
  expect_true(result_early$cum_hazard >= 0)
  expect_true(result_late$cum_hazard >= 0)

  # Late event should have higher cumulative hazard
  expect_true(result_late$cum_hazard > result_early$cum_hazard)
})

test_that(".solve_joint_ode biomarker values match observation times", {
  test_data <- list(
    time = 4,
    status = 1,
    longitudinal = list(
      times = c(0.5, 1.0, 2.0, 3.0),
      measurements = c(1.0, 1.2, 1.5, 1.7),
      n_obs = 4,
      covariates = NULL
    ),
    covariates = NULL
  )

  test_params <- list(
    coef = list(
      baseline = rep(0.1, 5),
      hazard = c(0.1, 0.05, 0.02),
      index_g = rep(0.05, 5),
      index_beta = c(0.01, 0.02, 0.005)  # Small values
    ),
    config = list(
      baseline = list(
        degree = 2,
        knots = c(1, 2),
        boundary_knots = c(0, 5),
        df = 5
      ),
      index = list(
        degree = 2,
        knots = c(0, 1),
        boundary_knots = c(-1, 2),
        df = 5
      )
    )
  )

  result <- JointODE:::.solve_joint_ode(test_data, test_params)

  # Biomarker values should be returned for each observation time
  expect_length(result$biomarker, test_data$longitudinal$n_obs)

  # Values should be finite and reasonable
  expect_true(all(is.finite(result$biomarker)))
})
