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
#' @param control List of optimization control parameters:
#'   \describe{
#'     \item{\code{method}}{Optimization algorithm (default: "L-BFGS-B")}
#'     \item{\code{maxit}}{Maximum iterations per M-step (default: 1000)}
#'     \item{\code{em_maxit}}{Maximum EM iterations (default: 10)}
#'     \item{\code{em_tol}}{EM convergence tolerance (default: 1e-4)}
#'     \item{\code{tol}}{Convergence tolerance (default: 1e-2)}
#'     \item{\code{verbose}}{Print progress and parameter estimates
#'       (default: FALSE)}
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
#' fit <- JointODE(
#'   longitudinal_formula = v ~ x1 + x2,
#'   longitudinal_data = sim$data$longitudinal_data,
#'   survival_formula = Surv(time, status) ~ w1 + w2,
#'   survival_data = sim$data$survival_data
#' )
#' summary(fit)
#' }
#'
#' @export
JointODE <- function(
  longitudinal_formula,
  longitudinal_data,
  survival_formula,
  survival_data,
  id = "id",
  time = "time",
  spline_baseline = list(
    degree = 3,
    n_knots = 5,
    knot_placement = "quantile",
    boundary_knots = NULL
  ),
  control = list(),
  ...
) {
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
    spline_baseline = spline_baseline
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
    function(s) s$longitudinal$n_obs > 0,
    data_process
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

  # 2.2 Initialize coefficients
  baseline_spline_coefficients <- rep(0, spline_baseline_config$df)
  hazard_coefficients <- rep(0, n_survival_covariates + 3)

  # Initialize beta without constraints (linear model)
  index_coefficients <- rnorm(n_longitudinal_covariates + 3, sd = 0.1)

  measurement_error_sd <- 1
  random_effect_sd <- 1

  # Set control defaults
  control_settings <- list(
    method = "L-BFGS-B",
    em_maxit = 100,
    maxit = 1000,
    tol = 1e-2,
    verbose = FALSE
  )
  control_settings[names(control)] <- control

  # EM Algorithm settings
  em_maxit <- control_settings$em_maxit %||% 10
  em_tol <- control_settings$em_tol %||% 1e-4

  # Build initial parameter structure
  parameters <- list(
    coefficients = list(
      baseline = baseline_spline_coefficients,
      hazard = hazard_coefficients,
      acceleration = index_coefficients,
      measurement_error_sd = measurement_error_sd,
      random_effect_sd = random_effect_sd
    ),
    configurations = list(
      baseline = spline_baseline_config
    )
  )

  # EM Algorithm iterations
  converged <- FALSE
  old_loglik <- -Inf
  loglik_history <- numeric(em_maxit)

  # Pre-compute dimensions
  n_baseline <- spline_baseline_config$df
  n_hazard <- n_survival_covariates + 3

  if (control_settings$verbose) {
    cat("\nStarting EM algorithm (max iterations:", em_maxit, ")\n")
    cat(sprintf(
      "%-5s %-15s %-10s %-10s %-10s\n",
      "Iter",
      "Log-Lik",
      "Change",
      "Sigma_e",
      "Sigma_b"
    ))
    cat(rep("-", 60), "\n", sep = "")
  }

  for (em_iter in seq_len(em_maxit)) {
    # E-step: Compute posteriors given current parameters
    posteriors <- .compute_posteriors(data_process, parameters)

    # Create fixed parameters structure for M-step
    fixed_parameters <- list(
      measurement_error_sd = parameters$coefficients$measurement_error_sd,
      random_effect_sd = parameters$coefficients$random_effect_sd
    )

    # Define objective function for M-step
    objective_fn <- function(par) {
      .compute_objective_joint(
        params = par,
        data_list = data_process,
        posteriors = posteriors,
        configurations = list(baseline = spline_baseline_config),
        fixed_parameters = fixed_parameters
      )
    }

    # Define gradient function for M-step
    gradient_fn <- function(par) {
      .compute_gradient_joint(
        params = par,
        data_list = data_process,
        posteriors = posteriors,
        configurations = list(baseline = spline_baseline_config),
        fixed_parameters = fixed_parameters
      )
    }

    # Initial parameter vector for this iteration
    par_init <- c(
      parameters$coefficients$baseline,
      parameters$coefficients$hazard,
      parameters$coefficients$acceleration
    )

    # M-step: Optimize parameters given posteriors
    res <- optim(
      par = par_init,
      fn = objective_fn,
      gr = gradient_fn,
      method = control_settings$method,
      control = list(
        maxit = control_settings$maxit,
        trace = if (control_settings$verbose) 1 else 0,
        REPORT = if (control_settings$verbose) 1 else 100
      )
    )

    # Extract and update parameters efficiently
    parameters$coefficients$baseline <- res$par[1:n_baseline]
    parameters$coefficients$hazard <- res$par[
      (n_baseline + 1):(n_baseline + n_hazard)
    ]
    parameters$coefficients$acceleration <- res$par[
      (n_baseline + n_hazard + 1):length(res$par)
    ]

    # Update variance components
    posteriors_updated <- .compute_posteriors(data_process, parameters)
    sds <- .compute_sds(data_process, parameters, posteriors_updated)
    parameters$coefficients$measurement_error_sd <- sds$measurement_error_sd
    parameters$coefficients$random_effect_sd <- sds$random_effect_sd

    # Track convergence
    new_loglik <- -res$value
    loglik_change <- new_loglik - old_loglik
    loglik_history[em_iter] <- new_loglik

    if (control_settings$verbose) {
      cat(sprintf(
        "%-5d %-15.6f %-10.6f %-10.4f %-10.4f\n",
        em_iter,
        new_loglik,
        loglik_change,
        parameters$coefficients$measurement_error_sd,
        parameters$coefficients$random_effect_sd
      ))

      # Always print parameter estimates when verbose
      cat("  \u03b7:", sprintf("%.3f", parameters$coefficients$baseline), "\n")
      cat(
        "  \u03b1:",
        sprintf("%.3f", parameters$coefficients$hazard[1:3]),
        "\n"
      )
      if (n_hazard > 3) {
        cat(
          "  \u03c6:",
          sprintf("%.3f", parameters$coefficients$hazard[4:n_hazard]),
          "\n"
        )
      }
      cat(
        "  \u03b2:",
        sprintf("%.3f", parameters$coefficients$acceleration),
        "\n\n"
      )
    }

    # Check convergence
    if (abs(loglik_change) < em_tol && em_iter > 1) {
      converged <- TRUE
      if (control_settings$verbose) {
        cat(rep("-", 60), "\n", sep = "")
        cat("EM converged successfully after", em_iter, "iterations\n")
        cat("Final log-likelihood:", sprintf("%.6f\n", new_loglik))
      }
      break
    }

    # Check for non-improvement
    if (em_iter > 3 && loglik_change < 0 && control_settings$verbose) {
      cat("\nWarning: Log-likelihood decreased. Possible convergence issues.\n")
    }

    old_loglik <- new_loglik
  }

  if (!converged && control_settings$verbose) {
    cat(rep("-", 60), "\n", sep = "")
    cat("EM did not converge within", em_maxit, "iterations\n")
    cat("Final log-likelihood:", sprintf("%.6f\n", new_loglik))
    cat("Consider increasing em_maxit or adjusting initial values\n")
  }

  # Final posteriors with converged parameters
  posteriors <- .compute_posteriors(data_process, parameters)
  measurement_error_sd <- parameters$coefficients$measurement_error_sd
  random_effect_sd <- parameters$coefficients$random_effect_sd

  # Calculate log-likelihood
  log_lik <- -res$value
  n_params <- length(res$par) + 2 # +2 for variance components
  aic <- -2 * log_lik + 2 * n_params
  bic <- -2 * log_lik + n_params * log(n_subjects)

  # Calculate residuals
  residuals_long <- numeric(n_observations)
  idx <- 1
  for (i in seq_along(data_process)) {
    subject <- data_process[[i]]
    if (subject$longitudinal$n_obs > 0) {
      # Solve ODE for this subject
      ode_solution <- .solve_joint_ode(subject, parameters)
      # Calculate residuals
      for (j in seq_len(subject$longitudinal$n_obs)) {
        residuals_long[idx] <- subject$longitudinal$measurements[j] -
          ode_solution$biomarker[j] -
          posteriors$b[i]
        idx <- idx + 1
      }
    }
  }

  # Return fitted model
  structure(
    list(
      coefficients = list(
        baseline = parameters$coefficients$baseline,
        hazard = parameters$coefficients$hazard,
        acceleration = parameters$coefficients$acceleration,
        measurement_error_sd = measurement_error_sd,
        random_effect_sd = random_effect_sd
      ),
      spline_config = list(
        baseline = spline_baseline_config
      ),
      logLik = log_lik,
      AIC = aic,
      BIC = bic,
      convergence = list(
        converged = converged,
        em_iterations = if (converged) em_iter else em_maxit,
        message = if (converged) {
          paste("EM algorithm converged after", em_iter, "iterations")
        } else {
          paste("EM algorithm did not converge within", em_maxit, "iterations")
        }
      ),
      fitted = list(
        longitudinal = NULL,
        survival = NULL
      ),
      residuals = list(
        longitudinal = residuals_long,
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
