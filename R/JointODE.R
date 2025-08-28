#' Fit Joint Model for Longitudinal and Survival Data Using ODEs
#'
#' @description
#' Fits a joint model that simultaneously analyzes longitudinal biomarker
#' trajectories and time-to-event outcomes using ordinary differential
#' equations (ODEs). This approach captures complex temporal dynamics in
#' biomarker evolution while quantifying their association with survival
#' through shared parameters.
#'
#' @param longitudinal_formula Formula for the longitudinal submodel.
#'   Left side specifies the response; right side includes time and
#'   covariates (e.g., \code{v ~ x1}).
#' @param longitudinal_data Data frame containing longitudinal measurements.
#'   Must include multiple observations per subject with columns for
#'   subject ID, measurement times, and response values.
#' @param survival_formula Formula for the survival submodel.
#'   Must use \code{Surv(time, status)} on the left side;
#'   right side specifies baseline covariates.
#' @param survival_data Data frame containing survival/event data.
#'   Must have exactly one row per subject with event time and status.
#' @param id Character string naming the subject ID variable.
#'   Must exist in both datasets (default: \code{"id"}).
#' @param time Character string naming the time variable in
#'   longitudinal data (default: \code{"time"}).
#' @param spline_baseline List of B-spline configuration for baseline hazard:
#'   \describe{
#'     \item{\code{degree}}{Degree of B-spline basis (default: 3)}
#'     \item{\code{n_knots}}{Number of interior knots (default: 5)}
#'     \item{\code{knot_placement}}{Method for knot placement: "quantile"
#'       (based on event times) or "equal" (default: "quantile")}
#'     \item{\code{boundary_knots}}{Boundary knots, if NULL uses event range
#'       (default: NULL)}
#'   }
#' @param spline_index List of B-spline configuration for single index model:
#'   \describe{
#'     \item{\code{degree}}{Degree of B-spline basis (default: 3)}
#'     \item{\code{n_knots}}{Number of interior knots (default: 4)}
#'     \item{\code{knot_placement}}{Method for knot placement: "quantile"
#'       (based on index values) or "equal" (default: "quantile")}
#'     \item{\code{boundary_knots}}{Boundary knots, if NULL uses index range
#'       (default: NULL)}
#'   }
#' @param control List of optimization control parameters:
#'   \describe{
#'     \item{\code{method}}{Optimization algorithm (default: "BFGS")}
#'     \item{\code{maxit}}{Maximum iterations (default: 1000)}
#'     \item{\code{tol}}{Convergence tolerance (default: 1e-6)}
#'     \item{\code{verbose}}{Print progress (default: FALSE)}
#'   }
#' @param ... Additional arguments passed to fitting functions.
#'
#' @return Object of class \code{"JointODE"} containing:
#'   \describe{
#'     \item{\code{coefficients}}{Estimated model parameters including
#'       longitudinal, survival, and association parameters}
#'     \item{\code{logLik}}{Log-likelihood at convergence}
#'     \item{\code{AIC}}{Akaike Information Criterion}
#'     \item{\code{BIC}}{Bayesian Information Criterion}
#'     \item{\code{convergence}}{Optimization convergence details}
#'     \item{\code{fitted}}{Fitted values for both submodels}
#'     \item{\code{residuals}}{Model residuals}
#'     \item{\code{data}}{Original input data}
#'     \item{\code{call}}{Matched function call}
#'   }
#'
#' @details
#' The joint ODE model links longitudinal and survival processes through
#' shared parameters. The longitudinal trajectory is modeled using ODEs
#' to capture non-linear dynamics, while the survival hazard incorporates
#' features of the trajectory (level, slope, or cumulative burden).
#'
#' Model estimation uses maximum likelihood with numerical integration
#' over random effects via adaptive Gauss-Hermite quadrature.
#'
#' @note
#' Input data are automatically validated and processed before fitting.
#' For large datasets or complex ODE systems, consider adjusting control
#' parameters to improve convergence.
#'
#' @concept modeling
#'
#' @importFrom stats optim
#' @importFrom survival Surv
#'
#' @examples
#' \dontrun{
#' sim <- simulate()
#' fit <- JointODE(
#'   longitudinal_formula = v ~ x1 + x2,
#'   longitudinal_data = sim$longitudinal_data,
#'   survival_formula = Surv(time, status) ~ w1 + w2,
#'   survival_data = sim$survival_data
#' )
#' summary(fit)
#' }
#'
#' @export
JointODE <- function(
    longitudinal_formula, longitudinal_data, survival_formula, survival_data,
    id = "id", time = "time",
    spline_baseline = list(
      degree = 3,
      n_knots = 5,
      knot_placement = "quantile",
      boundary_knots = NULL
    ),
    spline_index = list(
      degree = 3,
      n_knots = 4,
      knot_placement = "quantile",
      boundary_knots = NULL
    ),
    control = list(), ...) {
  # Store call
  cl <- match.call()

  # 1. Preprocessing
  # 1.1 Validate inputs
  .validate(
    longitudinal_formula = longitudinal_formula,
    longitudinal_data = longitudinal_data,
    survival_formula = survival_formula,
    survival_data = survival_data,
    id = id,
    time = time,
    spline_baseline = spline_baseline,
    spline_index = spline_index
  )

  # 1.2 Process data
  data_process <- .process(
    longitudinal_formula = longitudinal_formula,
    longitudinal_data = longitudinal_data,
    survival_formula = survival_formula,
    survival_data = survival_data,
    id = id,
    time = time
  )

  # 1.3 Extract properties
  event_times <- vapply(data_process, `[[`, numeric(1), "time")
  n_subjects <- attr(data_process, "n_subjects")
  n_observations <- sum(sapply(data_process, function(s) s$longitudinal$n_obs))
  subjects_with_long <- Filter(
    function(s) s$longitudinal$n_obs > 0, data_process
  )
  n_longitudinal_covariates <- if (length(subjects_with_long) > 0) {
    ncol(subjects_with_long[[1]]$longitudinal$covariates)
  } else {
    0
  }
  n_survival_covariates <- ncol(data_process[[1]]$covariates)

  # 2. Initialize parameters
  # 2.1 Configure splines
  spline_baseline_config <- .get_spline_config(
    x = event_times,
    degree = spline_baseline$degree,
    n_knots = spline_baseline$n_knots,
    knot_placement = spline_baseline$knot_placement,
    boundary_knots = spline_baseline$boundary_knots
  )
  spline_baseline_config$boundary_knots[1] <- 0

  scores <- seq(-0.3, 0.3, length.out = 100)
  spline_index_config <- .get_spline_config(
    x = scores,
    degree = spline_index$degree,
    n_knots = spline_index$n_knots,
    knot_placement = spline_index$knot_placement,
    boundary_knots = spline_index$boundary_knots
  )

  # 2.2 Initialize coefficients
  baseline_spline_coefficients <- rep(0, spline_baseline_config$df)
  hazard_coefficients <- rep(0, n_survival_covariates + 3)
  index_spline_coefficients <- rep(0, spline_index_config$df)

  # Initialize beta as unit vector, then convert to spherical coordinates
  index_coefficients_init <- rnorm(n_longitudinal_covariates + 3)
  norm_factor <- sqrt(sum(index_coefficients_init^2))
  if (norm_factor > 1e-10) {
    index_coefficients_init <- index_coefficients_init / norm_factor
  } else {
    index_coefficients_init[1] <- 1 # Default to first unit vector
  }

  # Convert to spherical coordinates for optimization
  index_theta <- .beta_to_spherical(index_coefficients_init)
  index_coefficients <- index_coefficients_init # Keep beta for initial E-step

  measurement_error_sd <- 1
  random_effect_sd <- 1


  # Set control defaults
  control_settings <- list(
    method = "L-BFGS-B",
    maxit = 1000,
    tol = 1e-2,
    verbose = FALSE
  )
  control_settings[names(control)] <- control
  b_hat <- numeric(n_subjects)

  # Track convergence
  old_params <- NULL
  converged <- FALSE
  actual_iter <- 0

  # EM algorithm implementation

  for (iter in seq_len(control_settings$maxit)) {
    actual_iter <- iter
    if (control_settings$verbose) {
      print(sprintf("Iteration %d/%d", iter, control_settings$maxit))
    }

    # E-step
    # Use current beta (already converted from theta if from M-step)
    parameters <- list(
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

    residuals <- list()
    posteriors <- list()

    # Pre-allocate vectors for posterior quantities
    posterior_mean_b <- numeric(n_subjects)
    posterior_variance_b <- numeric(n_subjects)
    posterior_expectation_exp_b <- numeric(n_subjects)

    for (i in seq_len(n_subjects)) {
      subject_data <- data_process[[i]]
      subject_id <- names(data_process)[i]
      ode_solution <- .solve_joint_ode(subject_data, parameters)
      posteriors[[subject_id]] <- .compute_posterior_aghq(
        ode_solution = ode_solution,
        data = subject_data,
        b_hat_init = b_hat[i],
        measurement_error_sd = measurement_error_sd,
        random_effect_sd = random_effect_sd
      )
      if (subject_data$longitudinal$n_obs > 0) {
        residuals[[subject_id]] <- subject_data$longitudinal$measurements -
          ode_solution$biomarker
      } else {
        residuals[[subject_id]] <- numeric(0)
      }

      # Store posterior quantities in vectors
      posterior_mean_b[i] <- posteriors[[subject_id]]$b_hat
      posterior_variance_b[i] <- posteriors[[subject_id]]$v_hat
      posterior_expectation_exp_b[i] <- posteriors[[subject_id]]$exp_b
      b_hat[i] <- posterior_mean_b[i]
    }

    # Create posteriors structure for efficient access (all subjects)
    posteriors_summary <- list(
      mean_b = posterior_mean_b,
      variance_b = posterior_variance_b,
      expectation_exp_b = posterior_expectation_exp_b
    )

    # M-step

    res_index <- optim(
      par = index_theta,
      fn = .compute_objective_beta,
      grad = .compute_grad_beta_forward,
      data_list = data_process,
      posteriors = posteriors_summary,
      config = list(
        baseline = spline_baseline_config,
        index = spline_index_config
      ),
      fixed_params = list(
        baseline = baseline_spline_coefficients,
        hazard = hazard_coefficients,
        index_g = index_spline_coefficients,
        measurement_error_sd = measurement_error_sd
      ),
      method = "L-BFGS-B",
      control = list(
        maxit = 10,
        trace = if (control_settings$verbose) 1 else 0,
        REPORT = if (control_settings$verbose) 1 else 10
      )
    )
    index_theta <- res_index$par
    index_coefficients <- .spherical_to_beta(index_theta)

    # Validate beta constraint ||beta|| = 1
    beta_norm <- sqrt(sum(index_coefficients^2))
    if (abs(beta_norm - 1) > 1e-10) {
      warning(
        sprintf(
          "Beta constraint violation: ||beta|| = %.6f, expected 1.0", beta_norm
        )
      )
      index_coefficients <- index_coefficients / beta_norm  # Normalize
    }

    if (control_settings$verbose) {
      cat("Beta:", index_coefficients, "\n")
    }

    # Optimize hazard coefficients
    par <- c(
      baseline_spline_coefficients,
      hazard_coefficients,
      index_spline_coefficients
    )
    res_hazard <- optim(
      par = par,
      fn = .compute_objective_theta,
      gr = .compute_grad_theta_forward,
      data_list = data_process,
      posteriors = posteriors_summary,
      config = list(
        baseline = spline_baseline_config,
        index = spline_index_config
      ),
      fixed_params = list(
        index_g = index_spline_coefficients,
        index_beta = index_coefficients
      ),
      method = "L-BFGS-B",
      control = list(
        maxit = 10,
        trace = if (control_settings$verbose) 1 else 0,
        REPORT = if (control_settings$verbose) 1 else 10
      )
    )
    baseline_spline_coefficients <- res_hazard$par[
      1:spline_baseline_config$df
    ]
    hazard_coefficients <- res_hazard$par[
      (spline_baseline_config$df + 1):
        (spline_baseline_config$df + n_survival_covariates + 3)
    ]
    if (control_settings$verbose) {
      cat("Lambda0:", baseline_spline_coefficients, "\n")
      cat("alpha:", hazard_coefficients[1:3], "\n")
      cat("eta:", hazard_coefficients[-(1:3)], "\n")
    }

    # Optimize variance components
    measurement_error_variance <- 0
    random_effect_variance <- 0
    for (i in seq_len(n_subjects)) {
      subject_data <- data_process[[i]]
      subject_id <- names(data_process)[i]
      residuals_i <- residuals[[subject_id]]
      n_i <- length(residuals_i)
      if (n_i > 0) {
        measurement_error_variance <- measurement_error_variance +
          sum(residuals_i^2) + n_i * posteriors[[subject_id]]$v_hat
      }
      random_effect_variance <- random_effect_variance +
        posterior_mean_b[i]^2 + posterior_variance_b[i]
    }
    measurement_error_sd <- sqrt(
      measurement_error_variance / max(n_observations, 1)
    )
    random_effect_sd <- sqrt(random_effect_variance / n_subjects)
    if (control_settings$verbose) {
      cat(
        "Sigma_e:", measurement_error_sd,
        "Sigma_b:", random_effect_sd, "\n"
      )
    }

    # Check convergence
    if (!is.null(old_params)) {
      param_diff <- sqrt(
        sum((baseline_spline_coefficients - old_params$baseline)^2) +
          sum((hazard_coefficients - old_params$hazard)^2) +
          sum((index_coefficients - old_params$index_beta)^2) +
          (measurement_error_sd - old_params$sigma_e)^2 +
          (random_effect_sd - old_params$sigma_b)^2
      )

      if (!is.na(param_diff) && param_diff < control_settings$tol) {
        converged <- TRUE
        if (control_settings$verbose) {
          cat("Converged at iteration", iter, "\n")
        }
        break
      }
    }

    # Store current parameters for next iteration
    old_params <- list(
      baseline = baseline_spline_coefficients,
      hazard = hazard_coefficients,
      index_beta = index_coefficients,
      sigma_e = measurement_error_sd,
      sigma_b = random_effect_sd
    )
  }

  # Return fitted model
  structure(
    list(
      coefficients = list(
        baseline = baseline_spline_coefficients,
        hazard = hazard_coefficients,
        index_g = index_spline_coefficients,
        index_beta = index_coefficients,
        measurement_error_sd = measurement_error_sd,
        random_effect_sd = random_effect_sd
      ),
      spline_config = list(
        baseline = spline_baseline_config,
        index = spline_index_config
      ),
      logLik = NA_real_,
      AIC = NA_real_,
      BIC = NA_real_,
      convergence = list(
        converged = converged,
        iterations = actual_iter,
        message = if (converged) {
          "Algorithm converged"
        } else {
          "Maximum iterations reached"
        }
      ),
      fitted = list(
        longitudinal = NULL,
        survival = NULL
      ),
      residuals = list(
        longitudinal = residuals,
        martingale = NULL
      ),
      random_effects = list(
        estimates = vapply(posteriors, `[[`, numeric(1), "b_hat"),
        variances = vapply(posteriors, `[[`, numeric(1), "v_hat")
      ),
      data = data_process,
      control = control_settings,
      call = cl
    ),
    class = "JointODE"
  )
}
