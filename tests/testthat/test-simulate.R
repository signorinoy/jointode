test_that("simulate generates valid data structure", {
  sim <- JointODE::simulate(n_subjects = 10, seed = 123)

  # Check output structure
  expect_type(sim, "list")
  expect_named(sim, c("longitudinal_data", "survival_data", "state"))

  # Check data frames
  expect_s3_class(sim$longitudinal_data, "data.frame")
  expect_s3_class(sim$survival_data, "data.frame")
  expect_s3_class(sim$state, "data.frame")

  # Check required columns
  expect_true(all(
    c("id", "time", "observed", "biomarker", "velocity", "acceleration") %in%
      names(sim$longitudinal_data)
  ))
  expect_true(all(c("id", "time", "status", "b") %in% names(sim$survival_data)))
  expect_named(sim$state, c("biomarker", "velocity"))
})

test_that("simulate respects n_subjects parameter", {
  n <- 10
  sim <- JointODE::simulate(n_subjects = n, seed = 456)

  expect_equal(nrow(sim$survival_data), n)
  expect_equal(length(unique(sim$longitudinal_data$id)), n)
  expect_equal(nrow(sim$state), n)
})

test_that("simulate respects seed for reproducibility", {
  sim1 <- JointODE::simulate(n_subjects = 10, seed = 789)
  sim2 <- JointODE::simulate(n_subjects = 10, seed = 789)

  expect_identical(sim1$survival_data$time, sim2$survival_data$time)
  expect_identical(
    sim1$longitudinal_data$observed,
    sim2$longitudinal_data$observed
  )
})

test_that("simulate handles custom parameters correctly", {
  sim <- JointODE::simulate(
    n_subjects = 10,
    shared_sd = 1.2,
    longitudinal = list(
      value = -0.15,
      slope = -0.9,
      ref = 0.5,
      time = 0.1,
      covariates = c(-0.2, 0.1),
      initial = list(ref = 1.2, covariates = c(-0.4, 0.2), random_coef = 1.0),
      error_sd = 0.5
    ),
    survival = list(
      baseline = list(type = "weibull", shape = 2, scale = 100),
      value = 0.01,
      slope = 0.5,
      covariates = c(-0.5, 0.3)
    ),
    covariates = list(
      trt = list(type = "binary", prob = 0.6),
      age = list(type = "normal", mean = 0, sd = 1)
    ),
    maxt = 50,
    seed = 111
  )

  expect_equal(ncol(sim$longitudinal_data), 8)
  expect_true(all(sim$survival_data$time <= 50))
  expect_true(
    mean(sim$survival_data$trt) > 0.4 && mean(sim$survival_data$trt) < 0.8
  )
})

test_that("simulate validates input parameters", {
  # Test only key validation cases to speed up tests
  expect_error(JointODE::simulate(n_subjects = -5), "positive integer")
  expect_error(JointODE::simulate(n_subjects = 10, shared_sd = -1), "positive")
  expect_error(
    JointODE::simulate(n_subjects = 10, longitudinal = "not_a_list"),
    "list"
  )
  expect_error(
    JointODE::simulate(
      n_subjects = 10,
      survival = list(baseline = list(type = "unknown"))
    ),
    "weibull"
  )
  # Dimension mismatch
  expect_error(
    JointODE::simulate(
      n_subjects = 10,
      longitudinal = list(
        value = -0.2,
        slope = -0.6,
        ref = 140,
        time = 0,
        covariates = c(1, 2), # 2 covariates
        initial = list(ref = 170, covariates = c(1), random_coef = 12),
        error_sd = 4
      )
    ),
    "same length"
  )
})

test_that("simulate generates reasonable biomarker trajectories", {
  sim <- JointODE::simulate(n_subjects = 20, seed = 222)

  # Check biomarker values are finite
  expect_true(all(is.finite(sim$longitudinal_data$biomarker)))
  expect_true(all(is.finite(sim$longitudinal_data$observed)))

  # Check velocities and accelerations are computed
  expect_true(all(is.finite(sim$longitudinal_data$velocity)))
  expect_true(all(is.finite(sim$longitudinal_data$acceleration)))

  # Check measurement error adds noise
  residuals <- sim$longitudinal_data$observed - sim$longitudinal_data$biomarker
  expect_true(sd(residuals) > 0)
})

test_that("simulate generates valid survival times", {
  sim <- JointODE::simulate(n_subjects = 20, maxt = 5, seed = 333)

  # All times should be positive and within maxt
  expect_true(all(sim$survival_data$time > 0))
  expect_true(all(sim$survival_data$time <= 5))

  # Status should be 0 or 1
  expect_true(all(sim$survival_data$status %in% c(0, 1)))

  # Should have at least some events (censoring is not guaranteed)
  expect_true(any(sim$survival_data$status == 1))

  # Test with parameters that ensure censoring
  sim_cens <- JointODE::simulate(
    n_subjects = 20,
    maxt = 2,
    survival = list(
      baseline = list(type = "weibull", shape = 1.5, scale = 20),
      value = 0.1,
      slope = 0.2,
      covariates = c(0.1, -0.1, -0.1)
    ),
    seed = 333
  )
  # With these parameters, we should have some censoring
  expect_true(any(sim_cens$survival_data$status == 0))
})

test_that("simulate handles edge cases", {
  # Small sample size
  sim_small <- JointODE::simulate(n_subjects = 1, seed = 444)
  expect_equal(nrow(sim_small$survival_data), 1)

  # Large error variance
  sim_noisy <- JointODE::simulate(
    n_subjects = 10,
    longitudinal = list(
      value = -0.2,
      slope = -0.6,
      ref = 0,
      time = 0,
      covariates = numeric(0),
      initial = list(ref = 1.5, covariates = numeric(0), random_coef = 1.0),
      error_sd = 10
    ),
    covariates = list(),
    seed = 555
  )
  expect_true(sd(sim_noisy$longitudinal_data$observed) > 5)

  # No covariates
  sim_no_cov <- JointODE::simulate(
    n_subjects = 10,
    longitudinal = list(
      value = -0.2,
      slope = -0.6,
      ref = 0,
      time = 0,
      covariates = numeric(0),
      initial = list(ref = 1.5, covariates = numeric(0), random_coef = 1.0),
      error_sd = 0.4
    ),
    survival = list(
      baseline = list(type = "weibull", shape = 1.6, scale = 150),
      value = 0.08,
      slope = 0.4,
      covariates = numeric(0)
    ),
    covariates = list(),
    seed = 666
  )
  expect_equal(ncol(sim_no_cov$longitudinal_data), 6) # No covariate columns
})

test_that("simulate creates consistent longitudinal observations", {
  skip_on_cran() # Skip on CRAN due to loop checks
  sim <- JointODE::simulate(n_subjects = 10, seed = 777)

  # Each subject should have observations up to their event/censoring time
  for (id in unique(sim$longitudinal_data$id)) {
    long_times <- sim$longitudinal_data$time[sim$longitudinal_data$id == id]
    surv_time <- sim$survival_data$time[sim$survival_data$id == id]

    expect_true(max(long_times) <= surv_time + 1e-6)
    expect_true(min(long_times) >= 0)

    # Times should be ordered
    expect_equal(long_times, sort(long_times))
  }
})

test_that("simulate respects covariate distributions", {
  skip_on_cran() # Skip on CRAN due to large sample size
  n <- 100 # Further reduced for faster tests
  sim <- JointODE::simulate(
    n_subjects = n,
    longitudinal = list(
      value = -0.2,
      slope = -0.6,
      ref = 0,
      time = 0,
      covariates = c(-0.2, 0.1), # Must match number of covariates
      initial = list(ref = 1.5, covariates = c(-0.5, 0.3), random_coef = 1.0),
      error_sd = 0.4
    ),
    survival = list(
      baseline = list(type = "weibull", shape = 1.6, scale = 150),
      value = 0.08,
      slope = 0.4,
      covariates = c(-0.7, 0.3) # Must match number of covariates
    ),
    covariates = list(
      binary_var = list(type = "binary", prob = 0.3),
      normal_var = list(type = "normal", mean = 0, sd = 1)
    ),
    seed = 888
  )

  # Binary variable
  prop <- mean(sim$survival_data$binary_var)
  expect_true(abs(prop - 0.3) < 0.15) # Wider tolerance for smaller sample

  # Normal variable
  expect_true(abs(mean(sim$survival_data$normal_var) - 0) < 0.4)
  expect_true(abs(sd(sim$survival_data$normal_var) - 1) < 0.4)
})

# Tests for .create_example_data
test_that(".create_example_data works correctly", {
  skip_on_cran() # Skip on CRAN - internal function test
  # Single call to create example data
  example <- JointODE:::.create_example_data(n_subjects = 20, seed = 123)

  # Test structure
  expect_type(example, "list")
  expect_named(example, c("data", "init"))
  expect_type(example$data, "list")
  expect_type(example$init, "list")
  expect_named(example$init, c("coefficients", "configurations"))

  # Test data generation
  expect_true(all(
    c("longitudinal_data", "survival_data", "state") %in%
      names(example$data)
  ))
  expect_equal(nrow(example$data$survival_data), 20)
  expect_equal(nrow(example$data$state), 20)

  # Test coefficients
  coef <- example$init$coefficients
  expect_true(all(
    c(
      "baseline",
      "acceleration",
      "hazard",
      "measurement_error_sd",
      "random_effect_sd"
    ) %in%
      names(coef)
  ))
  expect_type(coef$baseline, "double")
  expect_type(coef$acceleration, "double")
  expect_type(coef$hazard, "double")
  expect_equal(coef$measurement_error_sd, 0.1)
  expect_equal(coef$random_effect_sd, 0.1)

  # Test configurations
  config <- example$init$configurations
  expect_named(config, c("baseline", "autonomous"))
  expect_type(config$baseline, "list")
  expect_true(config$autonomous)
  expect_true(all(
    c("degree", "knots", "boundary_knots", "df") %in% names(config$baseline)
  ))

  # Test reproducibility
  example2 <- JointODE:::.create_example_data(n_subjects = 20, seed = 123)
  expect_identical(
    example$data$survival_data$time,
    example2$data$survival_data$time
  )
  expect_identical(example$init$coefficients, example2$init$coefficients)
})
