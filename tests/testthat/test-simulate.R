test_that("simulate() returns correct structure", {
  # Generate simulated data
  sim_data <- simulate(n = 20, seed = 123, verbose = FALSE)

  # Check return type
  expect_type(sim_data, "list")
  expect_named(sim_data, c("data.long", "data.surv"))

  # Check data.long structure
  expect_s3_class(sim_data$data.long, "data.frame")
  expect_named(sim_data$data.long, c("id", "time", "v", "x1", "x2"))
  expect_gt(nrow(sim_data$data.long), 0)

  # Check data.surv structure
  expect_s3_class(sim_data$data.surv, "data.frame")
  expect_named(sim_data$data.surv, c("id", "time", "status", "w1", "w2"))
  expect_equal(nrow(sim_data$data.surv), 20)

  # Check id consistency
  expect_equal(
    sort(unique(sim_data$data.long$id)),
    sort(unique(sim_data$data.surv$id))
  )
})

test_that("simulate() respects parameter inputs", {
  # Test with custom parameters
  sim_data <- simulate(
    n = 10,
    alpha = c(0.5, 0.2, -0.1),
    beta = c(-0.2, -0.4, 0.3, 0.15, 0.1),
    phi = c(0.3, -0.2),
    weibull_shape = 2,
    weibull_scale = 10,
    sigma_b = 0.3,
    sigma_e = 0.05,
    seed = 456,
    verbose = FALSE
  )

  # Check correct number of subjects
  expect_equal(length(unique(sim_data$data.surv$id)), 10)
  expect_equal(nrow(sim_data$data.surv), 10)

  # Check data types
  expect_type(sim_data$data.long$id, "integer")
  expect_type(sim_data$data.long$time, "double")
  expect_type(sim_data$data.long$v, "double")
  expect_type(sim_data$data.long$x1, "double")
  expect_type(sim_data$data.long$x2, "double")

  expect_type(sim_data$data.surv$id, "integer")
  expect_type(sim_data$data.surv$time, "double")
  expect_type(sim_data$data.surv$status, "integer")
  expect_type(sim_data$data.surv$w1, "double")
  expect_type(sim_data$data.surv$w2, "double")
})

test_that("simulate() produces valid survival data", {
  sim_data <- simulate(n = 50, seed = 789, verbose = FALSE)

  # Check survival times are positive
  expect_true(all(sim_data$data.surv$time > 0))

  # Check status is binary
  expect_true(all(sim_data$data.surv$status %in% c(0, 1)))

  # Check event rate is reasonable (between 20% and 80%)
  event_rate <- mean(sim_data$data.surv$status)
  expect_gt(event_rate, 0.2)
  expect_lt(event_rate, 0.8)

  # Check w2 is binary
  expect_true(all(sim_data$data.surv$w2 %in% c(0, 1)))

  # Check w1 is continuous
  expect_true(is.numeric(sim_data$data.surv$w1))
  expect_true(length(unique(sim_data$data.surv$w1)) > 1)
})

test_that("simulate() produces valid longitudinal data", {
  sim_data <- simulate(n = 30, seed = 321, verbose = FALSE)

  # Check multiple observations per subject
  obs_per_subject <- table(sim_data$data.long$id)
  expect_true(all(obs_per_subject > 1))

  # Check time values are reasonable
  expect_true(all(sim_data$data.long$time >= 0))
  expect_true(all(sim_data$data.long$time <= max(sim_data$data.surv$time)))

  # Check covariate values
  # x1 should be positive (exponential decay)
  expect_true(all(sim_data$data.long$x1 > 0))
  expect_true(all(sim_data$data.long$x1 <= 1))

  # x2 should be bounded (sine function)
  expect_true(all(abs(sim_data$data.long$x2) <= 0.2))

  # Check v values are numeric and finite
  expect_true(all(is.finite(sim_data$data.long$v)))
})

test_that("simulate() handles edge cases", {
  # Small sample size
  sim_small <- simulate(n = 2, seed = 111, verbose = FALSE)
  expect_equal(nrow(sim_small$data.surv), 2)
  expect_gt(nrow(sim_small$data.long), 2)

  # Large sample size
  sim_large <- simulate(n = 100, seed = 222, verbose = FALSE)
  expect_equal(nrow(sim_large$data.surv), 100)
  expect_gt(nrow(sim_large$data.long), 100)
})

test_that("simulate() reproducibility with seed", {
  # Same seed should produce identical results
  sim1 <- simulate(n = 25, seed = 999, verbose = FALSE)
  sim2 <- simulate(n = 25, seed = 999, verbose = FALSE)

  expect_identical(sim1$data.long, sim2$data.long)
  expect_identical(sim1$data.surv, sim2$data.surv)

  # Different seeds should produce different results
  sim3 <- simulate(n = 25, seed = 888, verbose = FALSE)
  expect_false(identical(sim1$data.long$v, sim3$data.long$v))
})

test_that("simulate() parameter validation", {
  # Test invalid n
  expect_error(
    simulate(n = 0, verbose = FALSE),
    "n must be a positive integer"
  )
  expect_error(
    simulate(n = -5, verbose = FALSE),
    "n must be a positive integer"
  )

  # Test invalid alpha length
  expect_error(
    simulate(n = 10, alpha = c(0.3, 0.1), verbose = FALSE),
    "alpha must be a numeric vector of length 3"
  )

  # Test invalid beta length
  expect_error(
    simulate(n = 10, beta = c(-0.3, -0.5), verbose = FALSE),
    "beta must be a numeric vector of length 5"
  )

  # Test invalid phi length
  expect_error(
    simulate(n = 10, phi = 0.2, verbose = FALSE),
    "phi must be a numeric vector of length 2"
  )

  # Test invalid sigma values
  expect_error(
    simulate(n = 10, sigma_b = -0.5, verbose = FALSE),
    "sigma_b must be positive"
  )
  expect_error(
    simulate(n = 10, sigma_e = 0, verbose = FALSE),
    "sigma_e must be positive"
  )

  # Test invalid weibull parameters
  expect_error(
    simulate(n = 10, weibull_shape = -1, verbose = FALSE),
    "weibull_shape must be positive"
  )
  expect_error(
    simulate(n = 10, weibull_scale = 0, verbose = FALSE),
    "weibull_scale must be positive"
  )
})

test_that("simulate() verbose output works", {
  # Test verbose = TRUE
  expect_message(
    simulate(n = 5, seed = 555, verbose = TRUE),
    "Generating ODE trajectories"
  )

  # Test verbose = FALSE (no output)
  expect_silent(
    simulate(n = 5, seed = 666, verbose = FALSE)
  )
})
