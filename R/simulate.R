#' Simulate Joint ODE Model Data
#'
#' @description
#' Generates synthetic data from a sophisticated joint modeling framework that
#' seamlessly integrates longitudinal biomarker trajectories with survival
#' outcomes through ordinary differential equations (ODEs). This simulation
#' engine produces realistic clinical trial datasets where complex biomarker
#' dynamics govern event hazards via shared random effects and trajectory
#' features, capturing the intricate interplay between disease progression
#' and time-to-event processes.
#'
#' @param n Integer. Number of subjects to simulate. Larger cohorts provide
#'   more stable parameter estimates (default: 100).
#' @param alpha Numeric vector of length 2. Association parameters quantifying
#'   how trajectory features influence survival hazard:
#'   \code{[biomarker, velocity]}.
#'   Positive values indicate increased risk (default: c(0.6, 1.0)).
#' @param beta Numeric vector governing ODE dynamics (length 5). Controls
#'   biomarker trajectory evolution: \code{[biomarker, velocity,
#'   x1, x2, time]} (default: c(-1.0, -0.6, -0.8, 0.5, 0.4)).
#' @param phi Numeric vector of length 2. Baseline covariate effects modulating
#'   survival hazard independently of biomarker dynamics: \code{[w1, w2]}
#'   (default: c(0.8, -1.2)).
#' @param weibull_shape Numeric. Weibull shape parameter (\eqn{\kappa})
#'   characterizing baseline hazard evolution. Values > 1 yield increasing
#'   hazard (aging effect), < 1 decreasing hazard (selection effect), =
#'   1 constant hazard (exponential) (default: 1).
#' @param weibull_scale Numeric. Weibull scale parameter (\eqn{\theta})
#'   determining the characteristic event time. Larger values shift the
#'   hazard curve rightward (default: 8).
#' @param sigma_b Numeric. Standard deviation of subject-specific random effects
#'   capturing unobserved heterogeneity. Larger values increase between-subject
#'   variability in both trajectories and hazards (default: 0.1).
#' @param sigma_e Numeric. Measurement error standard deviation reflecting
#'   assay precision and biological fluctuations. Smaller values indicate
#'   more reliable biomarker measurements (default: 0.1).
#' @param seed Integer. Random seed ensuring reproducible simulations.
#'   Essential for method validation and comparison studies (default: 42).
#' @param verbose Logical. Display informative progress messages during
#'   simulation workflow (default: FALSE).
#'
#' @return A list containing two complementary datasets:
#'   \describe{
#'     \item{longitudinal_data}{Data frame with longitudinal measurements:
#'       \itemize{
#'         \item \code{id}: Subject identifier
#'         \item \code{time}: Measurement time
#'         \item \code{v}: Observed biomarker value incorporating latent
#'           trajectory, random effect, and measurement error
#'         \item \code{x1}: First longitudinal covariate (standardized)
#'         \item \code{x2}: Second longitudinal covariate (standardized)
#'         \item \code{biomarker}: True latent biomarker trajectory value
#'         \item \code{velocity}: True latent biomarker velocity value
#'         \item \code{acceleration}: True latent biomarker acceleration value
#'       }}
#'     \item{survival_data}{Data frame with survival outcomes:
#'       \itemize{
#'         \item \code{id}: Subject identifier
#'         \item \code{time}: Observed event/censoring time
#'         \item \code{status}: Event indicator (1 = event, 0 = censored)
#'         \item \code{w1}: First survival covariate (standardized)
#'         \item \code{w2}: Second survival covariate (standardized)
#'         \item \code{b}: Shared random effect
#'       }}
#'   }
#'
#' @details
#' \strong{Mathematical Framework}
#'
#' The simulation orchestrates a sophisticated joint model architecture
#' comprising two intricately coupled components:
#'
#' \strong{1. Longitudinal Dynamics}
#'
#' Biomarker evolution is governed by a second-order linear ODE system
#' that captures temporal dynamics:
#' \deqn{\ddot{m}_i(t) = \boldsymbol{\beta}^{\top} \mathbf{Z}_i(t)}
#'
#' where:
#' \itemize{
#'   \item \eqn{m_i(t)}: Latent biomarker trajectory for subject \eqn{i}
#'   \item \eqn{\mathbf{Z}_i(t) = [m_i(t), \dot{m}_i(t), X_1, X_2, t]^{\top}}:
#'     State vector (5-dimensional)
#'   \item \eqn{\boldsymbol{\beta}}: Parameter vector governing homeostatic
#'     feedback, damping forces, and external influences
#' }
#'
#' Observed measurements arise from a hierarchical structure incorporating
#' both systematic and stochastic components:
#' \deqn{V_{ij} = m_i(T_{ij}) + b_i + \varepsilon_{ij}}
#'
#' where \eqn{b_i \sim N(0, \sigma_b^2)} captures subject-specific deviations
#' and \eqn{\varepsilon_{ij} \sim N(0, \sigma_e^2)} represents measurement
#' variability.
#'
#' \strong{2. Survival Process}
#'
#' The instantaneous hazard function elegantly links biomarker dynamics
#' to event risk through a multiplicative model:
#' \deqn{\lambda_i(t|b_i) = \lambda_0(t) \exp[\boldsymbol{\alpha}^{\top}
#'   \mathbf{m}_i(t) + \mathbf{W}_i^{\top} \boldsymbol{\phi} + b_i]}
#'
#' where:
#' \itemize{
#'   \item \eqn{\lambda_0(t) = (\kappa/\theta)(t/\theta)^{\kappa-1}}:
#'     Weibull baseline hazard capturing population-level risk evolution
#'   \item \eqn{\mathbf{m}_i(t) = [m_i(t), \dot{m}_i(t)]^{\top}}:
#'     Comprehensive trajectory feature vector
#'   \item \eqn{\mathbf{W}_i}: Time-invariant baseline characteristics
#'   \item \eqn{b_i}: Shared random effect inducing correlation between
#'     longitudinal and survival processes
#' }
#'
#' \strong{Simulation Workflow}
#'
#' The data generation proceeds through a carefully orchestrated pipeline:
#'
#' \enumerate{
#'   \item \strong{ODE Integration}: Numerical solution of the coupled ODE
#'     system using adaptive Runge-Kutta methods with automatic step size
#'     control for optimal accuracy-efficiency tradeoff
#'   \item \strong{Event Time Generation}: Sophisticated sampling from the
#'     conditional hazard distribution via the \code{simsurv} engine,
#'     accounting for complex time-varying covariates
#'   \item \strong{Censoring Mechanism}: Realistic administrative censoring
#'     uniformly distributed between the 50th and 95th percentiles of
#'     event times, mimicking clinical trial follow-up patterns
#'   \item \strong{Visit Scheduling}: Adaptive measurement protocols with
#'     higher frequency during critical periods (quarterly initially,
#'     semi-annual subsequently) and stochastic 10% missingness
#' }
#'
#' \strong{Parameter Interpretation Guide}
#'
#' \itemize{
#'   \item \strong{Trajectory Dynamics (\code{beta})}:
#'     \itemize{
#'       \item \code{beta[1]}: Homeostatic feedback strength (negative values
#'         promote stability)
#'       \item \code{beta[2]}: Damping coefficient controlling oscillation
#'         suppression
#'       \item \code{beta[3-4]}: Sensitivity to longitudinal covariates
#'       \item \code{beta[5]}: Time trend capturing systematic changes
#'     }
#'   \item \strong{Hazard Association (\code{alpha})}:
#'     \itemize{
#'       \item \code{alpha[1]}: Current biomarker value effect (positive
#'         indicates deleterious biomarker)
#'       \item \code{alpha[2]}: Velocity effect capturing prognostic value
#'         of trajectory direction
#'     }
#' }
#'
#' @note
#' \itemize{
#'   \item Visit schedules intelligently adapt to individual follow-up
#'     durations, balancing information gain with practical constraints
#'   \item Progress reporting can be silenced via \code{verbose = FALSE}
#'     for batch simulations
#' }
#'
#' @examples
#' # Basic usage with default parameters
#' sim <- simulate()
#' str(sim)
#'
#' # Check data characteristics
#' cat("Event rate:", mean(sim$survival_data$status), "\n")
#' cat(
#'   "Observations per subject:",
#'   nrow(sim$longitudinal_data) / nrow(sim$survival_data), "\n"
#' )
#'
#' # Visualize trajectories with survival information
#' library(ggplot2)
#' library(dplyr)
#'
#' # Select subjects with different outcomes
#' plot_data <- sim$longitudinal_data %>%
#'   left_join(sim$survival_data[, c("id", "time", "status")],
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
#' @concept data-simulation
#'
#' @importFrom deSolve ode
#' @importFrom simsurv simsurv
#' @importFrom stats rnorm runif quantile
#'
#' @export
simulate <- function(
  n = 100,
  alpha = c(0.6, 1.0),
  beta = c(-1.0, -0.6, -0.8, 0.5, 0.4),
  phi = c(0.8, -1.2),
  weibull_shape = 1,
  weibull_scale = 8,
  sigma_b = 0.1,
  sigma_e = 0.1,
  seed = 42,
  verbose = FALSE
) {
  # Validate inputs
  if (!is.numeric(n) || n <= 0 || n != as.integer(n)) {
    stop("n must be a positive integer")
  }
  if (!is.numeric(alpha) || length(alpha) != 2) {
    stop("alpha must be a numeric vector of length 2")
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
  parameters <- list(
    beta = beta,
    alpha = alpha,
    phi = phi,
    sigma_b = sigma_b,
    sigma_e = sigma_e,
    weibull = c(shape = weibull_shape, scale = weibull_scale)
  )

  # Generate survival times
  if (verbose) {
    message("Step 1/2: Generating survival times...")
  }
  survival_data <- .generate_survival_data(n, parameters)

  # Generate longitudinal measurements
  if (verbose) {
    message("Step 2/2: Generating longitudinal data...")
  }
  longitudinal_data <- .generate_longitudinal_data(survival_data, parameters)

  cols <- c("id", "time", "status", "w1", "w2", "b")
  survival_data <- survival_data[, cols]

  # Format and return
  list(longitudinal_data = longitudinal_data, survival_data = survival_data)
}

#' Compute Biomarker Acceleration
#'
#' @description
#' Computes the instantaneous acceleration (second derivative) of the
#' biomarker trajectory through linear combination of the current
#' state and covariate configuration.
#'
#' @param time Numeric. Current time point in the trajectory evolution.
#' @param biomarker Numeric. Instantaneous biomarker level \eqn{m(t)}.
#' @param velocity Numeric. Instantaneous rate of change \eqn{\dot{m}(t)}.
#' @param covariates Numeric vector. Subject-specific covariate values.
#' @param beta Numeric vector. Dynamics parameters.
#'
#' @return Numeric scalar. The computed acceleration \eqn{\ddot{m}(t)}
#'   representing the instantaneous rate of velocity change.
#'
#' @details
#' The acceleration is computed from the linear dynamics:
#' \deqn{\ddot{m}(t) = \boldsymbol{\beta}^{\top} \mathbf{z}(t)}
#' where the feature vector
#' \eqn{\mathbf{z}(t) = [m(t), \dot{m}(t), X_1, X_2, t]^{\top}}
#' encapsulates both endogenous state and exogenous influences.
#'
#' @noRd
.compute_acceleration_simulate <- function(
  time,
  biomarker,
  velocity,
  covariates,
  beta
) {
  z <- c(biomarker, velocity, covariates, time)
  sum(beta * z)
}

#' Solve Biomarker ODE System
#'
#' @description
#' Numerically integrates the coupled second-order ODE system to generate
#' complete biomarker trajectories including position, velocity, and
#' acceleration profiles over time.
#'
#' @param times Numeric vector. Evaluation time points for trajectory.
#' @param covariates Numeric vector. Subject-specific covariate profile.
#' @param parameters List containing model specifications:
#'   \itemize{
#'     \item \code{beta}: Dynamics parameters
#'   }
#'
#' @return Data frame containing the complete trajectory characterization:
#'   \itemize{
#'     \item \code{time}: Evaluation time points
#'     \item \code{biomarker}: Biomarker level \eqn{m(t)}
#'     \item \code{velocity}: Rate of change \eqn{\dot{m}(t)}
#'     \item \code{acceleration}: Rate of velocity change \eqn{\ddot{m}(t)}
#'   }
#'
#' @details
#' Employs the adaptive LSODA algorithm with automatic stiffness detection
#' and method switching for robust numerical integration. The system is
#' initialized at equilibrium: \eqn{m(0) = 0}, \eqn{\dot{m}(0) = 0}.
#'
#' @noRd
.solve_biomarker_ode <- function(times, covariates, parameters) {
  .biomarker_ode_deriv <- function(t, state, parms) {
    biomarker <- state[1]
    velocity <- state[2]
    acceleration <- .compute_acceleration_simulate(
      t,
      biomarker,
      velocity,
      covariates,
      parms$beta
    )
    list(c(velocity, acceleration))
  }
  ode_solution <- deSolve::ode(
    y = c(0, 0),
    times = sort(c(0, times)),
    func = .biomarker_ode_deriv,
    parms = parameters
  )
  idx <- match(times, ode_solution[, 1])
  biomarker <- ode_solution[idx, 2]
  velocity <- ode_solution[idx, 3]
  acceleration <- numeric(length(times))
  for (i in seq_along(times)) {
    acceleration[i] <- .compute_acceleration_simulate(
      times[i],
      biomarker[i],
      velocity[i],
      covariates,
      parameters$beta
    )
  }
  data.frame(
    time = times,
    biomarker = biomarker,
    velocity = velocity,
    acceleration = acceleration
  )
}

#' Generate Survival Data
#'
#' @description
#' Constructs a comprehensive survival dataset by simulating event times
#' from the joint model's complex hazard function, which seamlessly
#' integrates ODE-derived biomarker trajectories with baseline risks.
#'
#' @param n Integer. Target cohort size for simulation.
#' @param parameters List containing complete model specification:
#'   \itemize{
#'     \item \code{alpha}: Trajectory-hazard association parameters
#'     \item \code{phi}: Time-invariant covariate effects
#'     \item \code{weibull}: Baseline hazard shape and scale
#'     \item \code{sigma_b}: Random effect dispersion parameter
#'   }
#'
#' @return Data frame with columns:
#'   \itemize{
#'     \item \code{id}: Subject identifier
#'     \item \code{time}: Observed event/censoring time
#'     \item \code{status}: Event indicator (1 = event, 0 = censored)
#'     \item \code{w1}: First survival covariate (standardized)
#'     \item \code{w2}: Second survival covariate (standardized)
#'     \item \code{x1}: First longitudinal covariate (standardized)
#'     \item \code{x2}: Second longitudinal covariate (standardized)
#'     \item \code{b}: Subject-specific random effect
#'   }
#'
#' @details
#' Event generation leverages the sophisticated \code{simsurv} engine to
#' sample from hazard functions incorporating time-varying trajectory
#' features. Administrative censoring mimics realistic follow-up
#' limitations, uniformly distributed between empirical quantiles.
#'
#' @noRd
.generate_survival_data <- function(n, parameters) {
  # Extract baseline covariates
  covariates <- data.frame(
    id = seq_len(n),
    b = rnorm(n, 0, parameters$sigma_b),
    w1 = rnorm(n),
    w2 = rnorm(n),
    x1 = rnorm(n),
    x2 = rnorm(n)
  )

  # Define hazard function
  hazard_function <- function(t, x, betas, ...) {
    # Baseline hazard (Weibull)
    shape <- parameters$weibull["shape"]
    scale <- parameters$weibull["scale"]
    h0 <- (shape / scale) * (t / scale)^(shape - 1)
    biomarker <- .solve_biomarker_ode(
      t,
      as.numeric(c(x["x1"], x["x2"])),
      parameters
    )

    # Linear predictor
    eta <- parameters$alpha[1] *
      biomarker$biomarker +
      parameters$alpha[2] * biomarker$velocity +
      parameters$phi[1] * x["w1"] +
      parameters$phi[2] * x["w2"] +
      x["b"]
    h0 * exp(eta)
  }

  # Generate event times
  surv_times <- simsurv::simsurv(
    hazard = hazard_function,
    x = covariates,
    interval = c(0, 100)
  )

  # Apply censoring
  surv_times <- .apply_censoring(surv_times)

  # Combine with covariates
  merge(
    surv_times[, c("id", "time", "status")],
    covariates[, c("id", "w1", "w2", "x1", "x2", "b")],
    by = "id"
  )
}

#' Apply administrative censoring
#'
#' @description
#' Implements uniform administrative censoring to create
#' realistic censoring patterns in survival data.
#'
#' @param surv_times Data frame containing:
#'   \itemize{
#'     \item \code{eventtime}: True event times
#'     \item Other columns preserved from input
#'   }
#'
#' @return Data frame with censoring applied:
#'   \itemize{
#'     \item \code{time}: Minimum of event and censoring time
#'     \item \code{status}: 1 if event observed, 0 if censored
#'     \item All other input columns preserved
#'   }
#'
#' @details
#' Censoring times are uniformly distributed between the 50th
#' and 95th percentiles of event times, ensuring moderate
#' censoring rates typical of clinical trials.
#'
#' @noRd
.apply_censoring <- function(surv_times) {
  n <- nrow(surv_times)

  # Generate censoring times
  cens_times <- runif(
    n,
    min = quantile(surv_times$eventtime, 0.5),
    max = quantile(surv_times$eventtime, 0.95)
  )

  # Apply censoring
  surv_times$time <- pmin(surv_times$eventtime, cens_times)
  surv_times$status <- as.integer(surv_times$eventtime <= cens_times)
  surv_times
}

#' Generate longitudinal measurements
#'
#' @description
#' Creates longitudinal dataset with repeated measurements for all
#' subjects based on their follow-up times and ODE trajectories.
#'
#' @param survival_data Data frame with survival information including:
#'   \itemize{
#'     \item \code{id}: Subject identifiers
#'     \item \code{time}: Follow-up times
#'     \item \code{b}: Subject-specific random effects
#'   }
#' @param parameters List containing model parameters.
#'
#' @return Data frame with columns:
#'   \itemize{
#'     \item \code{id}: Subject identifier
#'     \item \code{time}: Measurement time
#'     \item \code{v}: Observed biomarker value
#'     \item \code{x1}: First longitudinal covariate (standardized)
#'     \item \code{x2}: Second longitudinal covariate (standardized)
#'     \item \code{biomarker}: True latent biomarker trajectory value
#'     \item \code{velocity}: True latent biomarker velocity value
#'     \item \code{acceleration}: True latent biomarker acceleration value
#'   }
#'
#' @details
#' Generates subject-specific visit schedules with 10% random
#' missingness to simulate realistic clinical trial patterns.
#' Measurements include ODE trajectory values plus random effects
#' and measurement error.
#'
#' @noRd
.generate_longitudinal_data <- function(survival_data, parameters) {
  n <- nrow(survival_data)
  data_list <- vector("list", n)
  for (i in seq_len(n)) {
    long_data_subject <- .generate_subject_measurements(
      obs_time = survival_data[i, "time"],
      covariates = c(
        survival_data[i, "x1"],
        survival_data[i, "x2"]
      ),
      random_effect = survival_data[i, "b"],
      parameters = parameters
    )
    long_data_subject$id <- survival_data[i, "id"]
    data_list[[i]] <- long_data_subject
  }
  long_data <- do.call(rbind, data_list)
  cols <- c(
    "id",
    "time",
    "v",
    "x1",
    "x2",
    "biomarker",
    "velocity",
    "acceleration"
  )
  long_data[, cols]
}

#' Generate measurements for one subject
#'
#' @description
#' Creates subject-specific longitudinal measurements following
#' a realistic visit schedule with random missingness.
#'
#' @param obs_time Numeric. Maximum observation time for the subject.
#' @param random_effect Numeric. Subject-specific random effect \eqn{b_i}.
#' @param parameters List containing model parameters.
#'
#' @return Data frame with measurements for one subject:
#'   \itemize{
#'     \item \code{time}: Visit times
#'     \item \code{v}: Observed biomarker values
#'     \item \code{x1}: First longitudinal covariate (standardized)
#'     \item \code{x2}: Second longitudinal covariate (standardized)
#'     \item \code{biomarker}: True latent biomarker trajectory value
#'     \item \code{velocity}: True latent biomarker velocity value
#'     \item \code{acceleration}: True latent biomarker acceleration value
#'   }
#'
#' @details
#' Visit schedule adapts to follow-up duration:
#' quarterly (0.25 years) for follow-up ≤ 2 years,
#' then quarterly for first 2 years plus semi-annual thereafter.
#' Applies 10% random missingness per visit.
#' Biomarker values: \eqn{V_{ij} = m_i(t_{ij}) + b_i + \varepsilon_{ij}}.
#'
#' @noRd
.generate_subject_measurements <- function(
  obs_time,
  covariates,
  random_effect,
  parameters
) {
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
    x_vals[i, ] <- covariates
  }

  biomarker <- .solve_biomarker_ode(visit_times, covariates, parameters)

  data.frame(
    time = visit_times,
    v = biomarker$biomarker +
      random_effect +
      rnorm(n_visits, 0, parameters$sigma_e),
    x1 = x_vals[, 1],
    x2 = x_vals[, 2],
    biomarker = biomarker$biomarker,
    velocity = biomarker$velocity,
    acceleration = biomarker$acceleration
  )
}

#' Create visit schedule based on observation time
#'
#' @description
#' Generates a realistic clinical visit schedule with varying
#' frequency based on follow-up duration.
#'
#' @param obs_time Numeric. Maximum observation time in years.
#'
#' @return Numeric vector of scheduled visit times.
#'
#' @details
#' Schedule follows clinical trial conventions:
#' \itemize{
#'   \item Quarterly visits (0.25 years) for short follow-up (≤ 2 years)
#'   \item Quarterly visits (every 0.25 years) for first 2 years when
#'           follow-up > 2 years
#'   \item Semi-annual visits (every 0.5 years) after 2 years
#' }
#' This pattern reflects intensive early monitoring followed by
#' maintenance phase visits.
#'
#' @noRd
.create_visit_schedule <- function(obs_time) {
  if (obs_time <= 2) {
    # Quarterly visits for short follow-up
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

#' Estimate B-spline Coefficients
#'
#' @description
#' Internal function that estimates B-spline coefficients to approximate
#' a given function using least squares fitting.
#'
#' @param x Numeric vector. Data points where the function is evaluated.
#' @param f0 Function. Target function to approximate with B-splines.
#' @param config List. B-spline configuration with components:
#'   \itemize{
#'     \item \code{degree}: Polynomial degree (default: 3)
#'     \item \code{n_knots}: Number of interior knots
#'     \item \code{knot_placement}: "quantile" or "equal"
#'     \item \code{boundary_knots}: Optional boundary knot locations
#'   }
#'
#' @return Numeric vector of B-spline coefficients that best approximate
#'   the target function at the given data points.
#'
#' @details
#' This function performs the following steps:
#' \enumerate{
#'   \item Configures the B-spline basis using the provided settings
#'   \item Evaluates the basis functions at data points
#'   \item Evaluates the target function at data points
#'   \item Solves the least squares problem to find optimal coefficients
#' }
#'
#' The resulting coefficients minimize
#' \eqn{||B\theta - f_0(x)||^2} where B is the basis matrix.
#'
#' @examples
#' \dontrun{
#' # Approximate a sine function
#' x <- seq(0, 2 * pi, length.out = 100)
#' f0 <- function(t) sin(t)
#' config <- list(degree = 3, n_knots = 5, knot_placement = "quantile")
#' coef <- .estimate_bspline_coef(x, f0, config)
#' }
#'
#' @noRd
.estimate_bspline_coef <- function(x, f0, config) {
  splin_config <- .get_spline_config(
    x = x,
    degree = config$degree,
    n_knots = config$n_knots,
    knot_placement = config$knot_placement,
    boundary_knots = config$boundary_knots
  )

  # Compute B-spline basis at grid points
  basis_matrix <- .compute_spline_basis(x, splin_config)

  # Compute target values
  y_target <- f0(x)

  # Ensure y_target is a column vector
  if (!is.matrix(y_target)) {
    y_target <- matrix(y_target, ncol = 1)
  }

  # Fit coefficients using least squares
  spline_coefficients <- solve(
    t(basis_matrix) %*% basis_matrix,
    t(basis_matrix) %*% y_target
  )
  as.vector(spline_coefficients)
}


#' Create Example Dataset for JointODE
#'
#' @description
#' Internal function that generates a complete example dataset for
#' demonstrating and testing the JointODE package. Creates simulated
#' data with known parameters for model validation.
#'
#' @param n Integer. Number of subjects to simulate (default: 100).
#'
#' @return A list containing:
#'   \describe{
#'     \item{data}{Simulated data from \code{simulate()} function:
#'       \itemize{
#'         \item \code{longitudinal_data}: Longitudinal measurements
#'         \item \code{survival_data}: Survival outcomes
#'       }
#'     }
#'     \item{formulas}{Model formulas:
#'       \itemize{
#'         \item \code{longitudinal}: v ~ x1 + x2
#'         \item \code{survival}: Surv(time, status) ~ w1 + w2
#'       }
#'     }
#'     \item{parameters}{True parameters and configurations:
#'       \itemize{
#'         \item \code{coefficients}: All model coefficients
#'         \item \code{configuration}: B-spline configurations
#'       }
#'     }
#'   }
#'
#' @details
#' This function serves as the data generation engine for the package's
#' example dataset. It:
#' \enumerate{
#'   \item Simulates data using predefined parameter values
#'   \item Estimates B-spline coefficients for baseline hazard and
#'     nonlinearity functions
#'   \item Packages everything into a structured format for analysis
#' }
#'
#' The generated dataset includes:
#' \itemize{
#'   \item Weibull baseline hazard (shape=1.5, scale=8)
#'   \item ODE dynamics with normalized parameters
#'   \item Measurement error SD = 0.1
#'   \item Random effect SD = 0.1
#'   \item B-spline approximations with degree 3
#' }
#'
#' @examples
#' \dontrun{
#' # Generate example data for 100 subjects (used for package's sim dataset)
#' example_data <- .create_example_data(n = 100)
#'
#' # Access components
#' head(example_data$data$longitudinal_data)
#' example_data$formulas$longitudinal
#' example_data$parameters$coefficients$hazard
#' }
#'
#' @seealso \code{\link{simulate}} for the underlying simulation engine
#'
#' @noRd
.create_example_data <- function(n = 100) {
  # Define default configurations (match JointODE defaults)
  spline_baseline <- list(
    degree = 3,
    n_knots = 5,
    knot_placement = "quantile",
    boundary_knots = NULL
  )

  # Define functions once (exponential baseline, Weibull shape=1)
  lambda_0 <- function(t) rep(log(1 / 8), length(t))

  # Generate and process data
  data <- simulate(n = n)

  # Extract times once
  times <- data$survival_data[, "time"]

  # Compute baseline spline coefficients
  baseline_spline_coefficients <- .estimate_bspline_coef(
    times,
    lambda_0,
    spline_baseline
  )

  # Define coefficients (use simulate function defaults)
  # hazard: [alpha1, alpha2, phi1, phi2]
  hazard_coefficients <- c(0.6, 1.0, 0.8, -1.2)
  # acceleration: [biomarker, velocity, intercept, x1, x2, time]
  # The intercept (3rd position) is added for estimation
  acceleration_coefficients <- c(-1.0, -0.6, 0, -0.8, 0.5, 0.4)

  # Create spline configurations
  spline_baseline_config <- .get_spline_config(
    x = times,
    degree = spline_baseline$degree,
    n_knots = spline_baseline$n_knots,
    knot_placement = spline_baseline$knot_placement,
    boundary_knots = spline_baseline$boundary_knots
  )

  # Return organized structure
  list(
    data = data,
    formulas = list(
      longitudinal = v ~ x1 + x2,
      survival = Surv(time, status) ~ w1 + w2
    ),
    parameters = list(
      coefficients = list(
        baseline = baseline_spline_coefficients,
        hazard = hazard_coefficients,
        acceleration = acceleration_coefficients,
        measurement_error_sd = 0.1,
        random_effect_sd = 0.1
      ),
      configurations = list(
        baseline = spline_baseline_config
      )
    )
  )
}
