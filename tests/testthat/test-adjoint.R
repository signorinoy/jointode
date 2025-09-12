test_that("adjoint requires jacobian_func", {
  ode_func <- function(t, x, params, data) {
    list(dx = -params[1] * x)
  }

  expect_error(
    adjoint(
      ode_func = ode_func,
      x0 = 1,
      params = 0.5,
      times = c(0, 1)
    ),
    "jacobian_func is required"
  )
})

test_that("adjoint requires objective_grad when objective_func provided", {
  ode_func <- function(t, x, params, data) {
    list(dx = -params[1] * x)
  }

  jacobian_func <- function(t, x, params, data) {
    list(
      df_dx = matrix(-params[1], 1, 1),
      df_dtheta = matrix(-x, 1, 1)
    )
  }

  objective_func <- function(x_final, data) {
    x_final^2
  }

  expect_error(
    adjoint(
      ode_func = ode_func,
      jacobian_func = jacobian_func,
      x0 = 1,
      params = 0.5,
      times = c(0, 1),
      objective_func = objective_func
    ),
    "objective_grad is required"
  )
})

test_that("adjoint computes correct gradients", {
  # Exponential decay system
  ode_func <- function(t, x, params, data) {
    list(dx = -params[1] * x)
  }

  jacobian_func <- function(t, x, params, data) {
    list(
      df_dx = matrix(-params[1], 1, 1),
      df_dtheta = matrix(-x, 1, 1)
    )
  }

  target <- 0.5
  objective_func <- function(x_final, data) {
    (x_final - target)^2
  }

  objective_grad <- function(x_final, data) {
    2 * (x_final - target)
  }

  # Test
  theta <- 0.5
  t_final <- 2

  result <- adjoint(
    ode_func = ode_func,
    jacobian_func = jacobian_func,
    x0 = 1,
    params = theta,
    times = c(0, t_final),
    objective_func = objective_func,
    objective_grad = objective_grad
  )

  # Analytical gradient
  x_t <- exp(-theta * t_final)
  dj_dtheta <- 2 * (x_t - target) * (-t_final) * x_t

  expect_equal(result$gradient, dj_dtheta, tolerance = 1e-7)
})

test_that("adjoint works with data parameter", {
  # ODE with external data
  ode_func <- function(t, x, params, data) {
    list(dx = -params[1] * data$A %*% x)
  }

  jacobian_func <- function(t, x, params, data) {
    list(
      df_dx = -params[1] * data$A,
      df_dtheta = matrix(-as.vector(data$A %*% x), ncol = 1)
    )
  }

  data <- list(A = matrix(c(1, 0.2, 0.2, 1), 2, 2))

  objective_func <- function(x_final, data) {
    sum(x_final^2)
  }

  objective_grad <- function(x_final, data) {
    2 * x_final
  }

  result <- adjoint(
    ode_func = ode_func,
    jacobian_func = jacobian_func,
    x0 = c(1, 1),
    params = 0.5,
    times = c(0, 2),
    data = data,
    objective_func = objective_func,
    objective_grad = objective_grad
  )

  expect_length(result$gradient, 1)
  expect_true(is.finite(result$gradient))
  expect_true(result$gradient < 0) # Increasing decay reduces final state
})

test_that("adjoint handles multi-parameter systems", {
  # 2D linear system with 2 parameters
  ode_func <- function(t, x, params, data) {
    dx1 <- -params[1] * x[1]
    dx2 <- -params[2] * x[2]
    list(dx = c(dx1, dx2))
  }

  jacobian_func <- function(t, x, params, data) {
    list(
      df_dx = diag(c(-params[1], -params[2])),
      df_dtheta = matrix(c(-x[1], 0, 0, -x[2]), 2, 2)
    )
  }

  objective_func <- function(x_final, data) {
    sum(x_final^2)
  }

  objective_grad <- function(x_final, data) {
    2 * x_final
  }

  # Test
  params <- c(0.3, 0.4)
  x0 <- c(1, 1.5)
  t_final <- 2

  result <- adjoint(
    ode_func = ode_func,
    jacobian_func = jacobian_func,
    x0 = x0,
    params = params,
    times = c(0, t_final),
    objective_func = objective_func,
    objective_grad = objective_grad
  )

  # Analytical gradients
  x1_t <- x0[1] * exp(-params[1] * t_final)
  x2_t <- x0[2] * exp(-params[2] * t_final)
  grad_analytical <- c(
    2 * x1_t * (-t_final) * x1_t,
    2 * x2_t * (-t_final) * x2_t
  )

  expect_equal(result$gradient, grad_analytical, tolerance = 1e-7)
})

test_that("adjoint with running cost", {
  # Simple system with running cost
  ode_func <- function(t, x, params, data) {
    list(dx = -params[1] * x)
  }

  jacobian_func <- function(t, x, params, data) {
    list(
      df_dx = matrix(-params[1], 1, 1),
      df_dtheta = matrix(-x, 1, 1)
    )
  }

  # Running cost: integral of x^2
  running_cost <- function(t, x, data) {
    x^2
  }

  running_cost_grad <- function(t, x, data) {
    2 * x
  }

  result <- adjoint(
    ode_func = ode_func,
    jacobian_func = jacobian_func,
    x0 = 1,
    params = 0.5,
    times = seq(0, 2, length.out = 21),
    running_cost = running_cost,
    running_cost_grad = running_cost_grad
  )

  # Check gradient is computed
  expect_true(is.finite(result$gradient))
  expect_true(result$gradient != 0)
})

test_that("adjoint validates Jacobian dimensions", {
  ode_func <- function(t, x, params, data) {
    list(dx = c(-params[1] * x[1], -params[2] * x[2]))
  }

  # Wrong dimensions for df_dx
  bad_jacobian <- function(t, x, params, data) {
    list(
      df_dx = matrix(-params[1], 1, 1), # Should be 2x2
      df_dtheta = matrix(c(-x[1], -x[2]), 2, 2)
    )
  }

  expect_error(
    adjoint(
      ode_func = ode_func,
      jacobian_func = bad_jacobian,
      x0 = c(1, 1),
      params = c(0.5, 0.3),
      times = c(0, 1)
    ),
    "df_dx must be"
  )
})

test_that("adjoint returns correct structure", {
  ode_func <- function(t, x, params, data) {
    list(dx = -params[1] * x)
  }

  jacobian_func <- function(t, x, params, data) {
    list(
      df_dx = matrix(-params[1], 1, 1),
      df_dtheta = matrix(-x, 1, 1)
    )
  }

  result <- adjoint(
    ode_func = ode_func,
    jacobian_func = jacobian_func,
    x0 = 1,
    params = 0.5,
    times = c(0, 1, 2),
    save_trajectory = TRUE
  )

  expect_s3_class(result, "adjoint")
  expect_named(
    result,
    c(
      "objective",
      "gradient",
      "final_state",
      "sensitivity_final",
      "trajectory",
      "n_states",
      "n_params",
      "times"
    )
  )
  expect_equal(result$n_states, 1)
  expect_equal(result$n_params, 1)
  expect_equal(nrow(result$trajectory), 3) # 3 time points
})
