#' Simulate Data from a Joint Ordinary Differential Equation Model
#'
#' @description
#' Generates synthetic longitudinal and time-to-event data under a joint
#' modeling
#' framework wherein the longitudinal biomarker trajectory is governed by a
#' second-order ordinary differential equation (ODE), and the survival process
#' is associated with the biomarker dynamics through shared random effects and
#' trajectory-dependent hazard functions. By default, biomarkers are generated
#' on a standardized scale (mean ~ 0, SD ~ 1) for better numerical stability.
#'
#' @param n_subjects Integer specifying the number of subjects to simulate
#'   (default: 500)
#' @param shared_sd Positive scalar defining the standard deviation of the
#'   shared
#'   random effects, \eqn{\sigma_b} (default: 0.1)
#' @param longitudinal List specifying the longitudinal sub-model parameters:
#'   \describe{
#'     \item{value}{Restoring force coefficient \eqn{\kappa} in the ODE
#'       system (default: -0.6)}
#'     \item{slope}{Damping coefficient \eqn{\gamma} controlling velocity
#'       feedback (default: -0.4)}
#'     \item{ref}{Reference (equilibrium) biomarker value \eqn{m_{ref}}
#'       (default: 0, for standardized biomarkers)}
#'     \item{time}{Time-varying coefficient for non-autonomous systems
#'       (default: 0, indicating autonomous dynamics)}
#'     \item{covariates}{Vector of fixed effects \eqn{\boldsymbol{\beta}}
#'       for covariates in the acceleration equation
#'       (default: c(0.8, -0.5, -0.5))}
#'     \item{initial}{List specifying initial condition parameters:
#'       \describe{
#'         \item{ref}{Baseline mean biomarker value \eqn{\mu_{ref}}
#'           (default: 1.0, for standardized biomarkers)}
#'         \item{covariates}{Vector of covariate effects
#'           \eqn{\boldsymbol{\beta}_{init}} on initial biomarker level
#'           (default: c(0.1, -0.1, 0.0))}
#'         \item{random_coef}{Scaling coefficient \eqn{\xi} for random effect
#'           influence on initial conditions (default: 0.0)}
#'       }
#'     }
#'     \item{error_sd}{Standard deviation \eqn{\sigma_{\epsilon}} of the
#'       measurement error process (default: 0.1, for standardized scale)}
#'     \item{n_measurements}{Number of longitudinal measurements per subject
#'       (default: 10)}
#'   }
#' @param survival List specifying the survival sub-model parameters:
#'   \describe{
#'     \item{baseline}{List defining the Weibull baseline hazard function:
#'       \describe{
#'         \item{type}{Character string specifying the baseline hazard type
#'           (currently only "weibull" is supported)}
#'         \item{shape}{Weibull shape parameter \eqn{k > 0}
#'           (default: 1.5)}
#'         \item{scale}{Weibull scale parameter \eqn{\lambda > 0}
#'           (default: 8.0)}
#'       }
#'     }
#'     \item{value}{Association parameter \eqn{\alpha_1} linking current
#'       biomarker value to hazard (default: 0.3)}
#'     \item{slope}{Association parameter \eqn{\alpha_2} linking biomarker
#'       velocity to hazard (default: 0.7)}
#'     \item{covariates}{Vector of regression coefficients
#'       \eqn{\boldsymbol{\phi}} for baseline covariates
#'       (default: c(0.4, -0.6))}
#'   }
#' @param covariates List defining the distributions of baseline covariates:
#'   \describe{
#'     \item{x1}{List with \code{type = "normal"}, \code{mean = 0},
#'       \code{sd = 1} for standardized continuous covariate}
#'     \item{x2}{List with \code{type = "normal"}, \code{mean = 0},
#'       \code{sd = 1} for standardized continuous covariate}
#'     \item{x3}{List with \code{type = "binary"} and \code{prob = 0.7}
#'       for binary covariate}
#'   }
#' @param maxt Positive scalar specifying the maximum follow-up time in the
#'   study (default: 10 days)
#' @param seed Integer seed for random number generation to ensure
#'   reproducibility (default: 42)
#'
#' @return A list containing three components:
#' \describe{
#'   \item{\code{longitudinal_data}}{Data frame comprising longitudinal
#'     observations with columns:
#'     \itemize{
#'       \item \code{id}: Subject identifier (integer)
#'       \item \code{time}: Observation time point (numeric)
#'       \item \code{observed}: Measured biomarker value including
#'         measurement error, \eqn{y_{ij}}
#'       \item \code{biomarker}: True underlying biomarker value,
#'         \eqn{m_i(t_{ij})}
#'       \item \code{velocity}: First derivative of the biomarker trajectory,
#'         \eqn{\dot{m}_i(t_{ij})}
#'       \item \code{acceleration}: Second derivative of the biomarker
#'         trajectory, \eqn{\ddot{m}_i(t_{ij})}
#'       \item \code{x1}, \code{x2}, \code{x3}:
#'         Time-invariant covariates
#'     }
#'   }
#'   \item{\code{survival_data}}{Data frame containing time-to-event data
#'     with columns:
#'     \itemize{
#'       \item \code{id}: Subject identifier
#'       \item \code{time}: Observed event or censoring time, \eqn{T_i}
#'       \item \code{status}: Event indicator, \eqn{\delta_i}
#'         (1 = event observed, 0 = censored)
#'       \item \code{b}: Realized subject-specific random effect, \eqn{b_i}
#'       \item \code{x1}, \code{x2}, \code{x3}:
#'         Baseline covariate values
#'     }
#'   }
#'   \item{\code{state}}{An \eqn{n \times 2} matrix containing initial states
#'     \eqn{[m_i(0), \dot{m}_i(0)]} for each subject}
#' }
#'
#' @details
#' The simulation framework implements a joint model comprising longitudinal
#' and survival
#' sub-models linked through shared random effects and trajectory-dependent
#' associations.
#'
#' \subsection{Longitudinal Sub-model}{
#' The biomarker trajectory \eqn{m_i(t)} for subject \eqn{i} evolves according
#' to the
#' second-order ordinary differential equation:
#' \deqn{\ddot{m}_i(t) = \kappa[m_i(t) - m_{ref}] + \gamma\dot{m}_i(t) +
#'   \mathbf{X}_i^T\boldsymbol{\beta} + \theta t}
#' where:
#' \itemize{
#'   \item \eqn{\kappa} represents the restoring force coefficient
#'   \item \eqn{\gamma} denotes the damping coefficient
#'   \item \eqn{m_{ref}} is the reference (equilibrium) value
#'   \item \eqn{\mathbf{X}_i} contains time-invariant covariates
#'   \item \eqn{\boldsymbol{\beta}} represents fixed effects
#'   \item \eqn{\theta} is the time-varying coefficient
#'     (0 for autonomous systems)
#' }
#'
#' Initial conditions are specified as:
#' \itemize{
#'   \item \eqn{m_i(0) = \mu_{ref} +
#'     \mathbf{X}_i^T\boldsymbol{\beta}_{init} + \xi b_i}
#'   \item \eqn{\dot{m}_i(0) = 0} (zero initial velocity)
#' }
#'
#' The observed longitudinal measurements incorporate additive Gaussian noise:
#' \deqn{y_{ij} = m_i(t_{ij}) + b_i + \epsilon_{ij}}
#' where \eqn{\epsilon_{ij} \sim \mathcal{N}(0, \sigma_\epsilon^2)}
#' represents independent measurement error.
#' }
#'
#' \subsection{Survival Sub-model}{
#' The instantaneous hazard function incorporates both current biomarker value
#' and velocity:
#' \deqn{\lambda_i(t) = \lambda_0(t)\exp(\alpha_1 m_i(t) + \alpha_2\dot{m}_i(t)
#'  + \mathbf{W}_i^T\boldsymbol{\phi} + b_i)}
#' where:
#' \itemize{
#'   \item \eqn{\lambda_0(t)} denotes the Weibull baseline hazard:
#'     \eqn{\lambda_0(t) = (k/\lambda)(t/\lambda)^{k-1}}
#'   \item \eqn{\alpha_1} quantifies the association with current biomarker
#'     value
#'   \item \eqn{\alpha_2} quantifies the association with biomarker velocity
#'   \item \eqn{\mathbf{W}_i} contains baseline covariates
#'   \item \eqn{\boldsymbol{\phi}} represents covariate effects on survival
#' }
#' }
#'
#' \subsection{Random Effects Structure}{
#' Subject-specific random effects \eqn{b_i \sim \mathcal{N}(0, \sigma_b^2)}
#' induce correlation
#' between the longitudinal and survival processes through:
#' \itemize{
#'   \item Initial biomarker level via \eqn{m_i(0)}
#'   \item Hazard function via the frailty term
#'   \item Longitudinal observations via additive shift
#' }
#' }
#'
#' @examples
#' # Generate a small dataset with default parameters
#' sim_data <- simulate(n_subjects = 50, seed = 123)
#'
#' # Examine the structure of generated data
#' str(sim_data)
#'
#' \dontrun{
#' # Customized simulation with modified longitudinal dynamics
#' sim_data <- simulate(
#'   n_subjects = 200,
#'   longitudinal = list(
#'     value = -2.5,      # Stronger restoring force
#'     slope = -4.0,      # Modified damping
#'     ref = -1,          # Lower equilibrium value (standardized)
#'     time = 0.1,        # Non-autonomous system
#'     covariates = c(-0.4, 0.2, 0.1),
#'     initial = list(
#'       ref = 2.5,       # Initial value (standardized)
#'       covariates = c(-1.0, 0.8, 0.4),
#'       random_coef = 1.0
#'     ),
#'     error_sd = 0.5     # Increased measurement error (standardized)
#'   )
#' )
#' }
#' @concept data-simulation
#'
#' @importFrom stats rnorm rbinom
#' @importFrom utils tail
#'
#' @export
simulate <- function(
  n_subjects = 200,
  shared_sd = 0.1,
  longitudinal = list(
    value = -0.6,
    slope = -0.4,
    ref = 0.0,
    time = 0.0,
    covariates = c(0.8, -0.5, -0.5),
    initial = list(
      ref = 1.0,
      covariates = c(0.1, -0.1, 0.0),
      random_coef = 0.0
    ),
    error_sd = 0.1,
    n_measurements = 10
  ),
  survival = list(
    baseline = list(
      type = "weibull",
      shape = 1.5,
      scale = 8.0
    ),
    value = 0.3,
    slope = 0.7,
    covariates = c(0.4, -0.6, -0.3)
  ),
  covariates = list(
    x1 = list(type = "normal", mean = 0, sd = 1),
    x2 = list(type = "normal", mean = 0, sd = 1),
    x3 = list(type = "binary", prob = 0.7)
  ),
  maxt = 10,
  seed = 42
) {
  # Validate basic parameters
  stopifnot(
    "n_subjects must be a positive integer" = is.numeric(n_subjects) &&
      length(n_subjects) == 1 &&
      n_subjects > 0 &&
      n_subjects == round(n_subjects),
    "shared_sd must be positive" = is.numeric(shared_sd) &&
      length(shared_sd) == 1 &&
      shared_sd > 0,
    "maxt must be positive" = is.numeric(maxt) && length(maxt) == 1 && maxt > 0
  )

  # Validate longitudinal structure first
  stopifnot(
    "longitudinal must be a list" = is.list(longitudinal)
  )

  # Set default for n_measurements if not provided
  if (is.null(longitudinal$n_measurements)) {
    longitudinal$n_measurements <- 10
  }

  # Validate remaining longitudinal parameters
  stopifnot(
    "longitudinal parameters must be numeric" = is.numeric(
      longitudinal$value
    ) &&
      is.numeric(longitudinal$slope) &&
      is.numeric(longitudinal$time),
    "longitudinal$ref must be numeric" = is.numeric(longitudinal$ref),
    "longitudinal$error_sd must be positive" = is.numeric(
      longitudinal$error_sd
    ) &&
      longitudinal$error_sd > 0,
    "longitudinal$n_measurements must be a positive integer" = is.numeric(
      longitudinal$n_measurements
    ) &&
      longitudinal$n_measurements > 0 &&
      longitudinal$n_measurements == round(longitudinal$n_measurements),
    "longitudinal$covariates must be numeric vector" = is.numeric(
      longitudinal$covariates
    ) &&
      is.vector(longitudinal$covariates),
    "longitudinal$initial must be a list" = is.list(longitudinal$initial),
    "initial parameters must be numeric" = is.numeric(
      longitudinal$initial$ref
    ) &&
      is.numeric(longitudinal$initial$random_coef),
    "initial$covariates must be numeric vector" = is.numeric(
      longitudinal$initial$covariates
    ) &&
      is.vector(longitudinal$initial$covariates)
  )

  # Validate survival structure
  stopifnot(
    "survival must be a list" = is.list(survival),
    "survival$baseline must be a list" = is.list(survival$baseline),
    "Only weibull baseline hazard is supported" = survival$baseline$type ==
      "weibull",
    "Weibull parameters must be positive" = is.numeric(
      survival$baseline$shape
    ) &&
      survival$baseline$shape > 0 &&
      is.numeric(survival$baseline$scale) &&
      survival$baseline$scale > 0,
    "survival coefficients must be numeric" = is.numeric(survival$value) &&
      is.numeric(survival$slope),
    "survival$covariates must be numeric vector" = is.numeric(
      survival$covariates
    ) &&
      is.vector(survival$covariates)
  )

  # Validate covariates
  stopifnot("covariates must be a list" = is.list(covariates))
  for (cov_name in names(covariates)) {
    cov <- covariates[[cov_name]]
    if (!is.list(cov)) {
      stop(paste("Covariate", cov_name, "must be a list"))
    }
    if (!cov$type %in% c("binary", "normal")) {
      stop(paste("Covariate", cov_name, "type must be binary or normal"))
    }
    if (cov$type == "binary") {
      if (!is.numeric(cov$prob) || cov$prob < 0 || cov$prob > 1) {
        stop(paste("Covariate", cov_name, "prob must be in [0,1]"))
      }
    } else {
      if (!is.numeric(cov$mean)) {
        stop(paste("Covariate", cov_name, "mean must be numeric"))
      }
      if (!is.numeric(cov$sd) || cov$sd <= 0) {
        stop(paste("Covariate", cov_name, "sd must be positive"))
      }
    }
  }

  # Validate dimension consistency
  n_cov <- length(covariates)
  if (n_cov > 0) {
    stopifnot(
      "All covariate vectors must have same length as covariates list" = length(
        longitudinal$covariates
      ) ==
        n_cov &&
        length(longitudinal$initial$covariates) == n_cov &&
        length(survival$covariates) == n_cov
    )
  }

  set.seed(seed)
  b <- .generate_shared_effects(n_subjects, shared_sd)
  x <- .generate_covariates(n_subjects, covariates)
  init <- .compute_initial_biomarker(
    x,
    b,
    longitudinal$initial
  )
  surv <- .generate_survival_data(
    x,
    b,
    init,
    longitudinal,
    survival,
    maxt
  )
  long <- .generate_longitudinal_data(
    x,
    b,
    init,
    surv,
    longitudinal
  )

  surv_data <- merge(surv, x, by = "id")
  surv_data <- merge(surv_data, b, by = "id")
  names(surv_data)[names(surv_data) == "eventtime"] <- "time"

  long_data <- merge(long, x, by = "id")
  col_order <- c(
    "id",
    "time",
    "observed",
    "biomarker",
    "velocity",
    "acceleration",
    names(covariates)
  )
  long_data <- long_data[, col_order]

  init_final <- init[, c("biomarker", "velocity")]

  list(
    longitudinal_data = long_data,
    survival_data = surv_data,
    state = init_final
  )
}

.create_example_data <- function(n_subjects = 200, seed = 123) {
  data <- simulate(n_subjects = n_subjects, seed = seed)
  coef_args <- formals(simulate)
  f0 <- function(t) {
    baseline_type <- eval(coef_args$survival$baseline$type)
    res <- switch(
      baseline_type,
      weibull = {
        shape <- eval(coef_args$survival$baseline$shape)
        scale <- eval(coef_args$survival$baseline$scale)
        (shape / scale) * (t / scale)^(shape - 1)
      },
      stop(paste(
        "Unsupported baseline hazard type:",
        baseline_type
      ))
    )
    log(res)
  }
  spline_config <- formals(JointODE)$spline_baseline
  spline_config <- .get_spline_config(
    x = data$survival_data$time,
    degree = spline_config$degree,
    n_knots = spline_config$n_knots,
    knot_placement = spline_config$knot_placement,
    boundary_knots = spline_config$boundary_knots
  )
  spline_config$boundary_knots[1] <- 0
  baseline_coef <- .estimate_bspline_coef(
    x = data$survival_data$time,
    f0 = f0,
    config = spline_config
  )
  acceleration_coef <- c(
    eval(coef_args$longitudinal$value),
    eval(coef_args$longitudinal$slope),
    -eval(coef_args$longitudinal$ref) *
      eval(coef_args$longitudinal$value),
    eval(coef_args$longitudinal$covariates),
    if (eval(coef_args$longitudinal$time) != 0) {
      eval(coef_args$longitudinal$time)
    } else {
      NULL
    }
  )
  hazard_coef <- c(
    eval(coef_args$survival$value),
    eval(coef_args$survival$slope),
    eval(coef_args$survival$covariates)
  )
  parameters <- list(
    coefficients = list(
      baseline = baseline_coef,
      acceleration = acceleration_coef,
      hazard = hazard_coef,
      measurement_error_sd = eval(coef_args$longitudinal$error_sd),
      random_effect_sd = eval(coef_args$shared_sd)
    ),
    configurations = list(
      baseline = spline_config,
      autonomous = eval(coef_args$longitudinal$time) == 0
    )
  )
  list(data = data, init = parameters)
}

.estimate_bspline_coef <- function(x, f0, config) {
  # Compute B-spline basis at grid points
  basis_matrix <- .compute_spline_basis(x, config)

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

.generate_shared_effects <- function(n_subjects, shared_sd) {
  data.frame(
    id = seq_len(n_subjects),
    b = rnorm(n_subjects, mean = 0, sd = shared_sd)
  )
}

.generate_covariates <- function(n_subjects, covariates) {
  covariate_data <- data.frame(id = seq_len(n_subjects))
  for (cov_name in names(covariates)) {
    cov_info <- covariates[[cov_name]]
    if (cov_info$type == "binary") {
      covariate_data[[cov_name]] <- rbinom(n_subjects, 1, cov_info$prob)
    } else if (cov_info$type == "normal") {
      covariate_data[[cov_name]] <- rnorm(
        n_subjects,
        cov_info$mean,
        cov_info$sd
      )
    } else {
      stop(paste("Unsupported covariate type:", cov_info$type))
    }
  }
  covariate_data
}

.compute_initial_biomarker <- function(x, b, initial) {
  value <- initial$ref +
    as.matrix(x[, names(x) != "id"]) %*% initial$covariates +
    initial$random_coef * b[, names(b) != "id"]
  slope <- rep(0, nrow(x))
  data.frame(id = x$id, biomarker = value, velocity = slope)
}

.solve_simple_biomarker_ode <- function(
  times = seq(0, 100, by = 0.1),
  init = c(1, 0),
  configurations = list()
) {
  .calculate_acceleration <- function(t, biomarker, velocity, parms) {
    omega <- 2 * pi / parms$period
    -2 *
      parms$xi *
      omega *
      velocity -
      (omega^2) * biomarker +
      parms$k * (omega^2) * parms$excitation(t)
  }
  .ode_deriv <- function(t, state, parms) {
    biomarker <- state[1]
    velocity <- state[2]
    acceleration <- .calculate_acceleration(
      t,
      biomarker,
      velocity,
      parms
    )
    list(c(velocity, acceleration))
  }
  default_config <- list(
    xi = -0.2,
    period = 20,
    k = 1.0,
    excitation = function(t) 0
  )
  ode_parms <- modifyList(default_config, configurations)
  ode_solution <- deSolve::ode(
    y = init,
    times = sort(c(0, times)),
    func = .ode_deriv,
    parms = ode_parms
  )
  acceleration <- .calculate_acceleration(
    ode_solution[, 1],
    ode_solution[, 2],
    ode_solution[, 3],
    ode_parms
  )
  idx <- match(times, ode_solution[, 1])
  data.frame(
    time = times,
    biomarker = ode_solution[idx, 2],
    velocity = ode_solution[idx, 3],
    acceleration = acceleration[idx]
  )
}


.solve_biomarker_ode <- function(times, x, init, longitudinal) {
  .biomarker_ode_deriv <- function(t, state, parms) {
    biomarker <- state[1]
    velocity <- state[2]
    acceleration <- parms$offset +
      parms$value * biomarker +
      parms$slope * velocity +
      parms$time * t
    list(c(velocity, acceleration))
  }

  offset <- -longitudinal$ref *
    longitudinal$value +
    sum(x * longitudinal$covariates)

  ode_parms <- list(
    offset = offset,
    value = longitudinal$value,
    slope = longitudinal$slope,
    time = longitudinal$time
  )

  ode_solution <- deSolve::ode(
    y = c(init$biomarker, init$velocity),
    times = sort(c(0, times)),
    func = .biomarker_ode_deriv,
    parms = ode_parms
  )
  idx <- match(times, ode_solution[, 1])
  biomarker <- ode_solution[idx, 2]
  velocity <- ode_solution[idx, 3]
  acceleration <- rep(offset, length(times)) +
    longitudinal$value * biomarker +
    longitudinal$slope * velocity +
    longitudinal$time * times
  data.frame(
    time = times,
    biomarker = biomarker,
    velocity = velocity,
    acceleration = acceleration
  )
}

.generate_survival_data <- function(
  x,
  b,
  init,
  longitudinal,
  survival,
  maxt
) {
  # Define hazard function
  hazard_function <- function(t, x, betas, ...) {
    # baseline hazard
    h0 <- switch(
      survival$baseline$type,
      weibull = {
        shape <- survival$baseline$shape
        scale <- survival$baseline$scale
        (shape / scale) * (t / scale)^(shape - 1)
      },
      stop(paste("Unsupported baseline hazard type:", survival$baseline$type))
    )
    x_cov <- x[!names(x) %in% c("id", "b", "biomarker", "velocity")]
    init <- x[c("biomarker", "velocity")]

    biomarker <- .solve_biomarker_ode(
      t,
      x_cov,
      init,
      longitudinal
    )
    linpred <- survival$value *
      biomarker$biomarker +
      survival$slope * biomarker$velocity +
      sum(x_cov * survival$covariates) +
      x$b

    h0 * exp(linpred)
  }

  covs <- merge(x, b, by = "id")
  covs <- merge(covs, init, by = "id")

  # Generate event times
  simsurv::simsurv(
    hazard = hazard_function,
    x = covs,
    maxt = maxt
  )
}

.generate_longitudinal_data <- function(x, b, init, surv, longitudinal) {
  n <- nrow(surv)
  data_list <- vector("list", n)
  maxt <- max(surv$eventtime)
  for (i in seq_len(n)) {
    times <- seq(
      0,
      surv[i, "eventtime"],
      by = maxt / longitudinal$n_measurements
    )
    if (tail(times, 1) == surv[i, "eventtime"]) {
      times <- times[-length(times)]
    }

    biomarkers <- .solve_biomarker_ode(
      times,
      x[x$id == surv[i, "id"], names(x) != "id"],
      init[init$id == surv[i, "id"], names(init) != "id"],
      longitudinal
    )
    biomarkers$id <- surv[i, "id"]
    random_effect <- b[b$id == surv[i, "id"], "b"]
    measurement_error <- rnorm(
      nrow(biomarkers),
      mean = 0,
      sd = longitudinal$error_sd
    )

    biomarkers$observed <- biomarkers$biomarker +
      random_effect +
      measurement_error

    data_list[[i]] <- biomarkers
  }
  long <- do.call(rbind, data_list)
  long
}
