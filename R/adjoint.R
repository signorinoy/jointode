#' Adjoint Sensitivity Analysis for ODE Systems
#'
#' @description
#' Computes gradients of scalar-valued objective functions with respect to ODE
#' parameters using the adjoint sensitivity method. This implementation requires
#' analytical derivatives for optimal performance and accuracy.
#'
#' @details
#' For a dynamical system described by ordinary differential equations:
#' \deqn{dx/dt = f(t, x, \theta, data)}
#' \deqn{x(t_0) = x_0}
#'
#' And an objective functional:
#' \deqn{J = g(x(T), data) + \int_{t_0}^{T} L(t, x(t), data) dt}
#'
#' The adjoint method efficiently computes the gradient \eqn{dJ/d\theta} by
#' solving:
#' 1. Forward ODE with sensitivity equations (forward pass)
#' 2. Adjoint ODE backward in time (if running cost present)
#'
#' This is particularly efficient when the number of parameters exceeds
#' the number of objective functions.
#'
#' @param ode_func Function defining the ODE system:
#'   function(t, x, params, data) returning list with element 'dx'
#'   containing the state derivatives dx/dt
#' @param jacobian_func Required function providing analytical derivatives:
#'   function(t, x, params, data) returning list with:
#'   - df_dx: State Jacobian \eqn{\partial f/\partial x} (n_states ×
#'     n_states matrix)
#'   - df_dtheta: Parameter Jacobian \eqn{\partial f/\partial \theta}
#'     (n_states × n_params matrix)
#' @param x0 Initial state vector at time t0
#' @param params Parameter vector \eqn{\theta} to compute sensitivities for
#' @param times Time grid for ODE integration (must include t0 and T)
#' @param data Optional list containing constant auxiliary data
#' @param objective_func Terminal cost function g(x(T), data)
#'   Returns scalar objective value at final time
#' @param objective_grad Required when objective_func is provided:
#'   function(x_final, data) returning gradient \eqn{\partial g/\partial x}
#'   at final time
#' @param running_cost Integrand function L(t, x, data) for running cost
#'   Returns scalar cost rate at time t
#' @param running_cost_grad Required when running_cost is provided:
#'   function(t, x, data) returning gradient \eqn{\partial L/\partial x}
#' @param rtol Relative error tolerance for ODE solver (default: 1e-8)
#' @param atol Absolute error tolerance for ODE solver (default: 1e-10)
#' @param method ODE solver algorithm (default: "lsoda" - adaptive solver)
#' @param save_trajectory Whether to return the full state trajectory
#'
#' @return Object of class "adjoint_result" containing:
#' \describe{
#'   \item{objective}{Scalar value of the objective function J}
#'   \item{gradient}{Gradient vector \eqn{dJ/d\theta} with respect to
#'     parameters}
#'   \item{final_state}{System state at final time x(T)}
#'   \item{sensitivity_final}{Sensitivity matrix
#'     \eqn{\partial x(T)/\partial \theta} at final time}
#'   \item{trajectory}{Full state trajectory (if save_trajectory = TRUE)}
#'   \item{n_states}{Number of state variables}
#'   \item{n_params}{Number of parameters}
#'   \item{times}{Time grid used for integration}
#' }
#'
#' @examples
#' \dontrun{
#' # Example: Parameter estimation for exponential decay
#'
#' # Define the ODE system: dx/dt = -\theta x
#' ode_system <- function(t, x, params, data) {
#'   list(dx = -params[1] * x)
#' }
#'
#' # Provide analytical Jacobians
#' jacobians <- function(t, x, params, data) {
#'   list(
#'     df_dx = matrix(-params[1], 1, 1),      # \partial f/\partial x
#'     df_dtheta = matrix(-x, 1, 1)           # \partial f/\partial \theta
#'   )
#' }
#'
#' # Define objective: squared error from target
#' target_value <- 0.5
#' objective <- function(x_final, data) {
#'   (x_final - target_value)^2
#' }
#'
#' # Gradient of objective
#' objective_gradient <- function(x_final, data) {
#'   2 * (x_final - target_value)
#' }
#'
#' # Compute sensitivity
#' result <- adjoint(
#'   ode_func = ode_system,
#'   jacobian_func = jacobians,
#'   x0 = 1,
#'   params = 0.5,
#'   times = seq(0, 2, length.out = 21),
#'   objective_func = objective,
#'   objective_grad = objective_gradient
#' )
#'
#' print(result)
#' }
#' @export
adjoint <- function(
  ode_func,
  jacobian_func,
  x0,
  params,
  times,
  data = NULL,
  objective_func = NULL,
  objective_grad = NULL,
  running_cost = NULL,
  running_cost_grad = NULL,
  rtol = 1e-8,
  atol = 1e-10,
  method = "lsoda",
  save_trajectory = FALSE
) {
  # ============================================================
  # Input Validation
  # ============================================================

  if (missing(jacobian_func) || is.null(jacobian_func)) {
    stop(
      "jacobian_func is required. ",
      "Please provide analytical Jacobians for optimal performance."
    )
  }

  if (!is.null(objective_func) && is.null(objective_grad)) {
    stop("objective_grad is required when objective_func is specified.")
  }

  if (!is.null(running_cost) && is.null(running_cost_grad)) {
    stop("running_cost_grad is required when running_cost is specified.")
  }

  if (!requireNamespace("deSolve", quietly = TRUE)) {
    stop(
      "Package 'deSolve' is required. ",
      "Please install it using: install.packages('deSolve')"
    )
  }

  # System dimensions
  n_states <- length(x0)
  n_params <- length(params)
  n_times <- length(times)

  # ============================================================
  # Jacobian Validation Helper
  # ============================================================

  validate_jacobian_matrices <- function(jac_list, n_states, n_params) {
    if (
      !is.list(jac_list) || !all(c("df_dx", "df_dtheta") %in% names(jac_list))
    ) {
      stop(
        "jacobian_func must return a list with ",
        "components 'df_dx' and 'df_dtheta'"
      )
    }

    if (
      !is.matrix(jac_list$df_dx) ||
        nrow(jac_list$df_dx) != n_states ||
        ncol(jac_list$df_dx) != n_states
    ) {
      stop(sprintf(
        "df_dx must be a %d \u00d7 %d matrix (got %d \u00d7 %d)",
        n_states,
        n_states,
        nrow(jac_list$df_dx),
        ncol(jac_list$df_dx)
      ))
    }

    if (
      !is.matrix(jac_list$df_dtheta) ||
        nrow(jac_list$df_dtheta) != n_states ||
        ncol(jac_list$df_dtheta) != n_params
    ) {
      stop(sprintf(
        "df_dtheta must be a %d \u00d7 %d matrix (got %d \u00d7 %d)",
        n_states,
        n_params,
        nrow(jac_list$df_dtheta),
        ncol(jac_list$df_dtheta)
      ))
    }
  }

  # ============================================================
  # Forward Pass: Augmented ODE System
  # ============================================================

  augmented_ode_rhs <- function(t, augmented_state, parms) {
    # Decompose augmented state vector
    state <- augmented_state[1:n_states]
    sensitivity_matrix <- matrix(
      augmented_state[(n_states + 1):(n_states + n_states * n_params)],
      nrow = n_states,
      ncol = n_params
    )

    # Evaluate ODE right-hand side
    ode_eval <- ode_func(t, state, parms, data)
    state_derivative <- if (is.list(ode_eval)) ode_eval$dx else ode_eval

    # Compute Jacobian matrices
    jacobian_matrices <- jacobian_func(t, state, parms, data)

    # Validate Jacobians (only on first evaluation for efficiency)
    if (t == times[1]) {
      validate_jacobian_matrices(jacobian_matrices, n_states, n_params)
    }

    # Forward sensitivity equation: dS/dt = (∂f/∂x)S + (∂f/∂θ)
    sensitivity_derivative <- jacobian_matrices$df_dx %*%
      sensitivity_matrix +
      jacobian_matrices$df_dtheta

    # Accumulate running cost if present
    cost_rate <- if (!is.null(running_cost)) {
      running_cost(t, state, data)
    } else {
      0
    }

    # Combine all derivatives
    list(c(
      as.vector(state_derivative),
      as.vector(sensitivity_derivative),
      cost_rate
    ))
  }

  # Set up initial conditions: [x0, S0 = 0, I0 = 0]
  augmented_initial <- c(
    x0, # Initial state
    rep(0, n_states * n_params), # Initial sensitivity (zero)
    0 # Initial cost integral
  )

  # Integrate forward ODE system
  forward_solution <- deSolve::ode(
    y = augmented_initial,
    times = times,
    func = augmented_ode_rhs,
    parms = params,
    rtol = rtol,
    atol = atol,
    method = method
  )

  # Extract components at final time
  final_time_idx <- n_times
  state_final <- as.vector(
    forward_solution[final_time_idx, 2:(n_states + 1)]
  )
  sensitivity_final <- matrix(
    forward_solution[
      final_time_idx,
      (n_states + 2):(n_states + 1 + n_states * n_params)
    ],
    nrow = n_states,
    ncol = n_params
  )
  cost_integral <- forward_solution[final_time_idx, ncol(forward_solution)]

  # Store trajectory if requested
  full_trajectory <- if (save_trajectory) {
    forward_solution[, 2:(n_states + 1), drop = FALSE]
  } else {
    NULL
  }

  # ============================================================
  # Objective Function Evaluation
  # ============================================================

  total_objective <- 0

  # Add terminal cost contribution
  if (!is.null(objective_func)) {
    total_objective <- total_objective + objective_func(state_final, data)
  }

  # Add integrated running cost
  if (!is.null(running_cost)) {
    total_objective <- total_objective + cost_integral
  }

  # ============================================================
  # Gradient Computation via Adjoint Method
  # ============================================================

  if (is.null(running_cost)) {
    # Case 1: Terminal cost only (direct computation)

    # Terminal adjoint: \u03bb(T) = \u2202g/\u2202x(T)
    adjoint_terminal <- if (!is.null(objective_grad)) {
      objective_grad(state_final, data)
    } else {
      rep(0, n_states)
    }

    # Gradient via chain rule: dJ/d\u03b8 = \u03bb(T)^T \u00b7 S(T)
    parameter_gradient <- as.vector(
      crossprod(adjoint_terminal, sensitivity_final)
    )
  } else {
    # Case 2: Including running cost (requires backward integration)

    # Store forward trajectory for interpolation
    state_trajectory <- forward_solution[, 2:(n_states + 1), drop = FALSE]

    # Define backward adjoint dynamics
    adjoint_dynamics <- function(backward_time, adjoint_vars, parms) {
      # Map backward time to forward time
      current_time <- parms$final_time - backward_time

      # Extract adjoint state
      # (parameter adjoint accumulated but not used directly)
      adjoint_state <- adjoint_vars[1:n_states]

      # Interpolate forward state at current time
      time_index <- findInterval(current_time, times)

      if (time_index == 0) {
        interpolated_state <- state_trajectory[1, ]
      } else if (time_index >= n_times) {
        interpolated_state <- state_trajectory[n_times, ]
      } else {
        # Linear interpolation weight
        alpha <- (current_time - times[time_index]) /
          (times[time_index + 1] - times[time_index])
        interpolated_state <- (1 - alpha) *
          state_trajectory[time_index, ] +
          alpha * state_trajectory[time_index + 1, ]
      }

      # Evaluate Jacobians at interpolated point
      jac_eval <- jacobian_func(current_time, interpolated_state, params, data)

      # Adjoint differential equations
      adjoint_state_deriv <- -crossprod(jac_eval$df_dx, adjoint_state)
      adjoint_param_deriv <- -crossprod(jac_eval$df_dtheta, adjoint_state)

      # Include running cost gradient contribution
      if (!is.null(running_cost_grad)) {
        running_gradient <- running_cost_grad(
          current_time,
          interpolated_state,
          data
        )
        adjoint_state_deriv <- adjoint_state_deriv - running_gradient
      }

      list(c(
        as.vector(adjoint_state_deriv),
        as.vector(adjoint_param_deriv)
      ))
    }

    # Set terminal conditions for adjoint system
    adjoint_terminal <- if (!is.null(objective_grad)) {
      objective_grad(state_final, data)
    } else {
      rep(0, n_states)
    }

    # Initial parameter gradient from terminal cost
    gradient_terminal <- as.vector(
      crossprod(adjoint_terminal, sensitivity_final)
    )

    # Set up backward integration interval
    time_final <- times[n_times]
    time_initial <- times[1]
    backward_interval <- c(0, time_final - time_initial)

    # Integrate adjoint system backward
    adjoint_solution <- deSolve::ode(
      y = c(adjoint_terminal, gradient_terminal),
      times = backward_interval,
      func = adjoint_dynamics,
      parms = list(final_time = time_final),
      rtol = rtol,
      atol = atol,
      method = method
    )

    # Extract gradient at initial time
    parameter_gradient <- adjoint_solution[
      2,
      (n_states + 2):(n_states + n_params + 1)
    ]
  }

  # ============================================================
  # Construct and Return Results
  # ============================================================

  structure(
    list(
      objective = as.numeric(total_objective),
      gradient = as.vector(parameter_gradient),
      final_state = state_final,
      sensitivity_final = sensitivity_final,
      trajectory = full_trajectory,
      n_states = n_states,
      n_params = n_params,
      times = times
    ),
    class = "adjoint_result"
  )
}

#' Print Method for Adjoint Results
#'
#' @param x An object of class "adjoint_result"
#' @param ... Additional arguments (unused)
#' @export
print.adjoint_result <- function(x, ...) {
  cat("Adjoint Sensitivity Analysis Results\n")
  cat("=====================================\n")
  cat(sprintf("Objective value: %.6g\n", x$objective))
  cat(sprintf(
    "Gradient: %s\n",
    paste(sprintf("%.6g", x$gradient), collapse = ", ")
  ))
  cat(sprintf(
    "Final state: %s\n",
    paste(sprintf("%.6g", x$final_state), collapse = ", ")
  ))
  cat(sprintf(
    "Problem dimensions: %d states, %d parameters\n",
    x$n_states,
    x$n_params
  ))
  invisible(x)
}
