# Test simulation functions when they are implemented in the package

test_that("parameter setup works correctly", {
  # Test parameter normalization
  beta <- c(-0.3, -0.5, 0.2, 0.1, 0.05)
  beta_norm <- beta / sqrt(sum(beta^2))
  
  expect_equal(sum(beta_norm^2), 1, tolerance = 1e-10)
  expect_true(all(is.finite(beta_norm)))
})

test_that("nonlinear link function behaves correctly", {
  # Test g(u) = 0.5 * tanh(u/3)
  g_function <- function(u) {
    0.5 * tanh(u / 3)
  }
  
  # Test bounds
  expect_lte(g_function(100), 0.5)
  expect_gte(g_function(-100), -0.5)
  
  # Test at zero
  expect_equal(g_function(0), 0)
  
  # Test monotonicity
  u_vals <- seq(-10, 10, by = 0.1)
  g_vals <- g_function(u_vals)
  expect_true(all(diff(g_vals) >= 0))  # Monotonically increasing
})

test_that("baseline hazard is valid", {
  baseline_hazard <- function(t, shape = 1.5, scale = 8) {
    (shape / scale) * (t / scale)^(shape - 1)
  }
  
  # Test at t = 0
  expect_equal(baseline_hazard(0, 1.5, 8), 0)
  
  # Test positivity for t > 0
  t_vals <- seq(0.1, 10, by = 0.1)
  h_vals <- baseline_hazard(t_vals, 1.5, 8)
  expect_true(all(h_vals > 0))
  
  # Test increasing hazard (shape > 1)
  expect_true(all(diff(h_vals) > 0))
})

test_that("random effects are generated correctly", {
  set.seed(123)
  n <- 100
  sigma_b <- 0.5
  
  b_values <- rnorm(n, 0, sigma_b)
  
  # Test distribution properties
  expect_equal(length(b_values), n)
  expect_equal(mean(b_values), 0, tolerance = 0.1)
  expect_equal(sd(b_values), sigma_b, tolerance = 0.1)
})

test_that("measurement error generation works", {
  set.seed(456)
  n_obs <- 1000
  sigma_e <- 0.1
  
  errors <- rnorm(n_obs, 0, sigma_e)
  
  # Test properties
  expect_equal(length(errors), n_obs)
  expect_equal(mean(errors), 0, tolerance = 0.01)
  expect_equal(sd(errors), sigma_e, tolerance = 0.02)
})

test_that("time-varying covariates are computed correctly", {
  # Treatment decay function
  x1_func <- function(t) exp(-t / 5)
  
  # Test at t = 0
  expect_equal(x1_func(0), 1)
  
  # Test decay
  expect_lt(x1_func(5), x1_func(0))
  expect_gt(x1_func(5), 0)
  
  # Seasonal effect function
  x2_func <- function(t) 0.2 * sin(2 * pi * t)
  
  # Test periodicity
  expect_equal(x2_func(0), x2_func(1), tolerance = 1e-10)
  expect_equal(x2_func(0.25), 0.2, tolerance = 1e-10)
  expect_equal(x2_func(0.75), -0.2, tolerance = 1e-10)
})

test_that("ODE initial conditions are valid", {
  m_0 <- 0
  m_dot_0 <- 0
  
  expect_true(is.finite(m_0))
  expect_true(is.finite(m_dot_0))
  expect_equal(length(m_0), 1)
  expect_equal(length(m_dot_0), 1)
})

test_that("data structure validation", {
  # Test longitudinal data structure
  test_long <- data.frame(
    id = c(1, 1, 2, 2),
    time = c(0, 1, 0, 1),
    value = c(0.5, 0.7, 0.3, 0.6),
    x1 = c(1, 0.8, 1, 0.8),
    x2 = c(0, 0.1, 0, 0.1),
    b_true = c(0.2, 0.2, -0.1, -0.1)
  )
  
  expect_true(is.data.frame(test_long))
  expect_true(all(c("id", "time", "value") %in% names(test_long)))
  expect_true(all(test_long$time >= 0))
  
  # Test survival data structure
  test_surv <- data.frame(
    id = c(1, 2),
    obstime = c(5.2, 3.8),
    status = c(1, 0),
    w1 = c(0.5, -0.3),
    w2 = c(1, 0)
  )
  
  expect_true(is.data.frame(test_surv))
  expect_true(all(c("id", "obstime", "status") %in% names(test_surv)))
  expect_true(all(test_surv$obstime > 0))
  expect_true(all(test_surv$status %in% c(0, 1)))
})