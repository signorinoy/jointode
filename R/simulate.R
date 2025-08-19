#' Simulate Joint ODE Model Data
#'
#' @description
#' Generates synthetic data from a joint ordinary differential equation (ODE)
#' model that couples longitudinal biomarker trajectories with time-to-event
#' outcomes. This framework produces realistic clinical trial or observational
#' study data where biomarker dynamics influence event hazards through shared
#' random effects and trajectory features.
#'
#' @param n Integer. Number of subjects to simulate (default: 50).
#' @param alpha Numeric vector of length 3. Association parameters linking
#'   trajectory features to survival hazard: \code{[value, velocity,
#'   acceleration]} (default: c(0.3, 0.1, -0.05)).
#' @param beta Numeric vector of length 5. ODE dynamics parameters controlling
#'   trajectory evolution: \code{[position_feedback, damping, treatment_effect,
#'   seasonal_effect, time_trend]} (default: c(-0.3, -0.5, 0.2, 0.1, 0.05)).
#' @param phi Numeric vector of length 2. Baseline covariate effects on
#'   survival hazard: \code{[continuous_covariate, binary_covariate]}
#'   (default: c(0.2, -0.15)).
#' @param weibull_shape Numeric. Shape parameter (\eqn{\kappa}) for Weibull
#'   baseline hazard. Values > 1 indicate increasing hazard, < 1 decreasing
#'   (default: 1.5).
#' @param weibull_scale Numeric. Scale parameter (\eqn{\theta}) for Weibull
#'   baseline hazard controlling the time scale (default: 8).
#' @param sigma_b Numeric. Standard deviation of subject-specific random
#'   intercepts representing between-subject heterogeneity (default: 0.5).
#' @param sigma_e Numeric. Standard deviation of measurement error for
#'   longitudinal observations (default: 0.1).
#' @param seed Integer. Random seed for reproducibility (default: 42).
#' @param verbose Logical. Print progress messages during simulation
#'   (default: TRUE).
#'
#' @return A list of containing:
#'   \describe{
#'     \item{data.long}{Data frame with longitudinal measurements:
#'       \itemize{
#'         \item \code{id}: Subject identifier
#'         \item \code{time}: Measurement time
#'         \item \code{v}: Observed biomarker value (latent + random effect
#'           + error)
#'         \item \code{x1}: Time-varying covariate 1 (treatment decay:
#'           \eqn{e^{-t/5}})
#'         \item \code{x2}: Time-varying covariate 2 (seasonal:
#'           \eqn{0.2\sin(2\pi t)})
#'       }}
#'     \item{data.surv}{Data frame with survival outcomes:
#'       \itemize{
#'         \item \code{id}: Subject identifier
#'         \item \code{time}: Observed event/censoring time
#'         \item \code{status}: Event indicator (1 = event, 0 = censored)
#'         \item \code{w1}: Baseline continuous covariate ~ N(0,1)
#'         \item \code{w2}: Baseline binary covariate ~ Bernoulli(0.5)
#'       }}
#'   }
#'
#' @details
#' \strong{Model Framework}
#'
#' The simulation implements a joint model with two coupled components:
#'
#' \strong{1. Longitudinal Model}
#'
#' Biomarker evolution follows a second-order nonlinear ODE:
#' \deqn{\ddot{m}_i(t) = g(\boldsymbol{\beta}^T \mathbf{Z}_i(t))}
#'
#' where:
#' \itemize{
#'   \item \eqn{m_i(t)} is the latent biomarker trajectory for subject \eqn{i}
#'   \item \eqn{\mathbf{Z}_i(t) = [m_i(t), \dot{m}_i(t), X_1(t), X_2(t), t]^T}
#'     is the augmented state vector
#'   \item \eqn{g(u) = 0.5 \tanh(u/3)} is a bounded nonlinear link function
#'     ensuring stability
#'   \item \eqn{\boldsymbol{\beta}} controls dynamics (feedback, damping,
#'     covariate effects)
#' }
#'
#' Observed measurements incorporate random effects and error:
#' \deqn{V_{ij} = m_i(T_{ij}) + b_i + \varepsilon_{ij}}
#'
#' with \eqn{b_i \sim N(0, \sigma_b^2)} and
#' \eqn{\varepsilon_{ij} \sim N(0, \sigma_e^2)}.
#'
#' \strong{2. Survival Model}
#'
#' The hazard function links trajectory dynamics to event risk:
#' \deqn{\lambda_i(t|b_i) = \lambda_0(t) \exp[\boldsymbol{\alpha}^T
#'   \mathbf{m}_i(t) + \mathbf{W}_i^T \boldsymbol{\phi} + b_i]}
#'
#' where:
#' \itemize{
#'   \item \eqn{\lambda_0(t) = (\kappa/\theta)(t/\theta)^{\kappa-1}} is the
#'     Weibull baseline hazard
#'   \item \eqn{\mathbf{m}_i(t) = [m_i(t), \dot{m}_i(t), \ddot{m}_i(t)]^T}
#'     captures trajectory features
#'   \item \eqn{\mathbf{W}_i} are baseline covariates
#'   \item \eqn{b_i} is the shared random effect linking models
#' }
#'
#' \strong{Simulation Process}
#'
#' \enumerate{
#'   \item \strong{ODE Integration}: Solve the second-order ODE system for each
#'     subject using adaptive LSODA algorithm
#'   \item \strong{Survival Generation}: Sample event times from the conditional
#'     hazard using the \code{simsurv} package
#'   \item \strong{Censoring}: Apply administrative censoring uniformly between
#'     50th and 95th percentiles of event times
#'   \item \strong{Longitudinal Sampling}: Generate visit schedules (quarterly
#'     for years 0-2, semi-annual thereafter) with 10% random missingness
#' }
#'
#' \strong{Parameter Interpretation}
#'
#' \itemize{
#'   \item \strong{ODE parameters (\code{beta})}:
#'     \itemize{
#'       \item \code{beta[1]}: Negative feedback (homeostasis), typically < 0
#'       \item \code{beta[2]}: Damping coefficient (stability), typically < 0
#'       \item \code{beta[3-5]}: Covariate and time effects
#'     }
#'   \item \strong{Association (\code{alpha})}:
#'     \itemize{
#'       \item \code{alpha[1]}: Current value association (typically > 0 for
#'         risk biomarkers)
#'       \item \code{alpha[2]}: Rate of change association
#'       \item \code{alpha[3]}: Acceleration association (often < 0 for
#'         stability)
#'     }
#' }
#'
#' @note
#' \itemize{
#'   \item The \code{beta} parameters are normalized internally to ensure
#'     model identifiability
#'   \item Visit schedules automatically adjust based on follow-up duration
#'   \item Setting \code{verbose = FALSE} suppresses all progress messages
#' }
#'
#' @examples
#' # Basic usage with default parameters
#' sim <- simulate()
#' str(sim)
#'
#' # Check data characteristics
#' cat("Event rate:", mean(sim$data.surv$status), "\n")
#' cat(
#'   "Observations per subject:",
#'   nrow(sim$data.long) / nrow(sim$data.surv), "\n"
#' )
#'
#' # Visualize trajectories with survival information
#' library(ggplot2)
#' library(dplyr)
#'
#' # Select subjects with different outcomes
#' plot_data <- sim$data.long %>%
#'   left_join(sim$data.surv[, c("id", "time", "status")],
#'     by = "id", suffix = c("", "_event")
#'   ) %>%
#'   filter(id %in% sample(unique(id), 9))
#'
#' ggplot(plot_data, aes(x = time, y = v)) +
#'   geom_line(aes(color = factor(status)), alpha = 0.7) +
#'   geom_point(aes(color = factor(status)), size = 1) +
#'   geom_vline(aes(xintercept = time_event),
#'     linetype = "dashed", alpha = 0.5
#'   ) +
#'   facet_wrap(~id, scales = "free_y", ncol = 3) +
#'   scale_color_manual(
#'     values = c("0" = "blue", "1" = "red"), labels = c("Censored", "Event")
#'   ) +
#'   theme_minimal() +
#'   labs(
#'     x = "Time (years)", y = "Biomarker Value", color = "Outcome",
#'     title = "Simulated Longitudinal Trajectories"
#'   )
#'
#' @seealso
#' \code{\link[simsurv]{simsurv}} for the underlying survival simulation
#' \code{\link[deSolve]{ode}} for ODE integration details
#'
#' @concept utilities
#'
#' @importFrom deSolve ode
#' @importFrom simsurv simsurv
#' @importFrom stats rnorm rbinom rexp runif approxfun quantile
#'
#' @export
simulate <- function(
    n = 50, alpha = c(0.3, 0.1, -0.05), beta = c(-0.3, -0.5, 0.2, 0.1, 0.05),
    phi = c(0.2, -0.15), weibull_shape = 1.5, weibull_scale = 8,
    sigma_b = 0.5, sigma_e = 0.1, seed = 42, verbose = TRUE) {
  # Validate inputs
  if (!is.numeric(n) || n <= 0 || n != as.integer(n)) {
    stop("n must be a positive integer")
  }
  if (!is.numeric(alpha) || length(alpha) != 3) {
    stop("alpha must be a numeric vector of length 3")
  }
  if (!is.numeric(beta) || length(beta) != 5) {
    stop("beta must be a numeric vector of length 5")
  }
  if (!is.numeric(phi) || length(phi) != 2) {
    stop("phi must be a numeric vector of length 2")
  }
  if (!is.numeric(weibull_shape) || weibull_shape <= 0) {
    stop("weibull_shape must be positive")
  }
  if (!is.numeric(weibull_scale) || weibull_scale <= 0) {
    stop("weibull_scale must be positive")
  }
  if (!is.numeric(sigma_b) || sigma_b <= 0) {
    stop("sigma_b must be positive")
  }
  if (!is.numeric(sigma_e) || sigma_e <= 0) {
    stop("sigma_e must be positive")
  }

  set.seed(seed)

  # Prepare parameters
  params <- list(
    beta = beta / sqrt(sum(beta^2)), # Normalize for identifiability
    alpha = alpha,
    phi = phi,
    sigma_b = sigma_b,
    sigma_e = sigma_e,
    weibull = c(shape = weibull_shape, scale = weibull_scale)
  )

  # Generate trajectories with progress tracking
  if (verbose) message("Step 1/3: Generating ODE trajectories...")
  trajectories <- vector("list", n)
  for (i in seq_len(n)) {
    trajectories[[i]] <- .generate_trajectory(i, params)
  }

  # Generate survival times
  if (verbose) message("Step 2/3: Generating survival times...")
  data.surv <- .generate_survival_data(trajectories, params, verbose)

  # Generate longitudinal measurements
  if (verbose) message("Step 3/3: Generating longitudinal data...")
  data.long <- .generate_longitudinal_data(
    trajectories, data.surv$obstime, params
  )

  # Format and return
  list(
    data.long = data.long,
    data.surv = data.surv[, c("id", "time", "status", "w1", "w2")]
  )
}

# Internal helper functions ----

#' Generate individual trajectory
#' @noRd
.generate_trajectory <- function(i, params) {
  # Subject-specific parameters
  b_i <- rnorm(1, 0, params$sigma_b)
  w_i <- c(rnorm(1), rbinom(1, 1, 0.5))

  # Create covariate function once
  x_func <- .create_covariate_function()
  # Solve ODE system
  sol <- .solve_ode(params$beta, x_func)

  # Create trajectory functions
  list(
    id = i,
    b_i = b_i,
    w_i = w_i,
    m_func = approxfun(sol$time, sol$m, rule = 2),
    m_dot_func = approxfun(sol$time, sol$m_dot, rule = 2),
    m_ddot_func = approxfun(sol$time, sol$m_ddot, rule = 2),
    x_func = x_func
  )
}

#' Solve ODE system for trajectory
#' @noRd
.solve_ode <- function(beta, x_func = NULL, tmax = 10, dt = 0.05) {
  # Create covariate function if not provided
  if (is.null(x_func)) {
    x_func <- .create_covariate_function()
  }
  # Define ODE system
  ode_system <- function(t, state, parms) {
    x_t <- parms$x_func(t)
    z <- c(state, x_t, t) # [m, m_dot, x1, x2, t]
    u <- sum(parms$beta * z)
    m_ddot <- 0.5 * tanh(u / 3) # Bounded nonlinearity
    list(c(state[2], m_ddot))
  }

  # Solve ODE
  times <- seq(0, tmax, by = dt)
  sol <- deSolve::ode(
    y = c(m = 0, m_dot = 0),
    times = times,
    func = ode_system,
    parms = list(beta = beta, x_func = x_func),
    method = "lsoda"
  )

  # Vectorized computation of acceleration
  n_times <- length(times)
  x_vals <- matrix(NA, n_times, 2)
  for (i in seq_len(n_times)) {
    x_vals[i, ] <- x_func(times[i])
  }
  z_mat <- cbind(sol[, 2:3], x_vals, sol[, 1])
  m_ddot <- 0.5 * tanh(as.vector(z_mat %*% beta) / 3)

  list(
    time = sol[, 1],
    m = sol[, 2],
    m_dot = sol[, 3],
    m_ddot = m_ddot
  )
}

#' Create time-varying covariate function
#' @noRd
.create_covariate_function <- function() {
  function(t) {
    c(
      exp(-t / 5), # Treatment decay
      0.2 * sin(2 * pi * t) # Seasonal effect
    )
  }
}

#' Generate survival data
#' @noRd
.generate_survival_data <- function(trajectories, params, verbose) {
  n <- length(trajectories)

  # Extract baseline covariates
  covdat <- data.frame(
    id = seq_len(n),
    b_i = vapply(trajectories, `[[`, numeric(1), "b_i"),
    w1 = vapply(trajectories, function(x) x$w_i[1], numeric(1)),
    w2 = vapply(trajectories, function(x) x$w_i[2], numeric(1))
  )

  # Define hazard function
  hazard_func <- function(t, x, betas, ...) {
    traj <- trajectories[[as.integer(x["id"])]]

    # Baseline hazard (Weibull)
    shape <- params$weibull["shape"]
    scale <- params$weibull["scale"]
    h0 <- (shape / scale) * (t / scale)^(shape - 1)

    # Linear predictor
    eta <- params$alpha[1] * traj$m_func(t) +
      params$alpha[2] * traj$m_dot_func(t) +
      params$alpha[3] * traj$m_ddot_func(t) +
      x["b_i"] +
      betas["w1"] * x["w1"] +
      betas["w2"] * x["w2"]

    h0 * exp(eta)
  }

  # Generate event times
  surv_times <- simsurv::simsurv(
    hazard = hazard_func,
    x = covdat,
    betas = c(w1 = params$phi[1], w2 = params$phi[2]),
    maxt = 10,
    interval = c(1e-8, 15)
  )

  # Apply censoring
  surv_times <- .apply_censoring(surv_times, verbose)

  # Combine with covariates
  surv_times$time <- surv_times$obstime
  merge(surv_times, covdat[, c("id", "w1", "w2")], by = "id")
}

#' Apply administrative censoring
#' @noRd
.apply_censoring <- function(surv_times, verbose) {
  n <- nrow(surv_times)

  # Generate censoring times
  cens_times <- runif(n,
    min = quantile(surv_times$eventtime, 0.5),
    max = quantile(surv_times$eventtime, 0.95)
  )

  # Apply censoring
  surv_times$obstime <- pmin(surv_times$eventtime, cens_times)
  surv_times$status <- as.integer(surv_times$eventtime <= cens_times)

  if (verbose) {
    censoring_rate <- 100 * (1 - mean(surv_times$status))
    message(sprintf("  Achieved censoring: %.1f%%", censoring_rate))
  }

  surv_times
}

#' Generate longitudinal measurements
#' @noRd
.generate_longitudinal_data <- function(trajectories, obs_times, params) {
  # Pre-allocate list for efficiency
  n <- length(trajectories)
  data_list <- vector("list", n)
  for (i in seq_len(n)) {
    data_list[[i]] <- .generate_subject_measurements(
      trajectories[[i]],
      obs_times[i],
      params
    )
  }

  do.call(rbind, data_list)
}

#' Generate measurements for one subject
#' @noRd
.generate_subject_measurements <- function(traj, obs_time, params) {
  # Define visit schedule
  visit_times <- .create_visit_schedule(obs_time)

  # Apply random missingness (10% dropout)
  n_visits <- length(visit_times)
  keep <- runif(n_visits) > 0.1
  visit_times <- visit_times[keep]
  if (length(visit_times) == 0) {
    visit_times <- 0
    n_visits <- 1
  } else {
    n_visits <- length(visit_times)
  }

  # Vectorized computation
  x_vals <- matrix(NA, n_visits, 2)
  for (i in seq_len(n_visits)) {
    x_vals[i, ] <- traj$x_func(visit_times[i])
  }

  data.frame(
    id = rep(traj$id, n_visits),
    time = visit_times,
    v = traj$m_func(visit_times) + traj$b_i +
      rnorm(n_visits, 0, params$sigma_e),
    x1 = x_vals[, 1],
    x2 = x_vals[, 2]
  )
}

#' Create visit schedule based on observation time
#' @noRd
.create_visit_schedule <- function(obs_time) {
  if (obs_time <= 2) {
    # Frequent visits for short follow-up
    seq(0, obs_time, by = 0.25)
  } else {
    # Quarterly for first 2 years, then semi-annual
    quarterly <- seq(0, 2, by = 0.25)
    if (obs_time > 2.5) {
      semiannual <- seq(2.5, obs_time, by = 0.5)
      c(quarterly, semiannual)
    } else {
      quarterly
    }
  }
}
