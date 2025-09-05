#' Joint Modeling of Longitudinal and Survival Data Using ODEs
#'
#' @description
#' Implements a unified framework for jointly modeling longitudinal biomarker
#' trajectories and time-to-event outcomes using ordinary differential
#' equations (ODEs). The model captures complex non-linear dynamics in
#' biomarker evolution while simultaneously quantifying their association with
#' survival risk through
#' shared random effects and flexible hazard specifications.
#'
#' @param longitudinal_formula A formula specifying the longitudinal submodel.
#'   The left-hand side defines the response variable, while the right-hand
#'   side specifies fixed effects including time-varying and baseline
#'   covariates
#'   (e.g., \code{biomarker ~ time + treatment + age}).
#' @param longitudinal_data A data frame containing repeated measurements with
#'   one row per observation. Required columns include subject identifier,
#'   measurement times, response values, and any covariates specified in the
#'   formula.
#' @param survival_formula A formula for the survival submodel using
#'   \code{Surv(time, status)} notation on the left-hand side. The right-hand
#'   side specifies baseline hazard covariates
#'   (e.g., \code{Surv(event_time, event) ~ treatment + age}).
#' @param survival_data A data frame with time-to-event information containing
#'   one row per subject. Must include event/censoring times, event indicators,
#'   and baseline covariates.
#' @param id Character string specifying the column name for subject
#'   identifiers. This variable must be present in both longitudinal and
#'   survival datasets
#'   to link observations (default: \code{"id"}).
#' @param time Character string specifying the column name for measurement
#'   times
#'   in the longitudinal dataset (default: \code{"time"}).
#' @param spline_baseline A list controlling the B-spline representation of the
#'   baseline hazard function with the following components:
#'   \describe{
#'     \item{\code{degree}}{Polynomial degree of the B-spline basis functions
#'       (default: 3, cubic splines)}
#'     \item{\code{n_knots}}{Number of interior knots for flexibility
#'       (default: 5, providing moderate flexibility)}
#'     \item{\code{knot_placement}}{Strategy for positioning knots:
#'       \code{"quantile"} places knots at quantiles of observed event times,
#'       \code{"equal"} uses equally-spaced knots (default: \code{"quantile"})}
#'     \item{\code{boundary_knots}}{A numeric vector of length 2 specifying
#'       the boundary knot locations. If \code{NULL}, automatically set to the
#'       range of observed event times (default: \code{NULL})}
#'   }
#' @param control A list of optimization and algorithmic settings:
#'   \describe{
#'     \item{\code{method}}{Optimization algorithm for parameter estimation.
#'       Options include \code{"L-BFGS-B"}, \code{"BFGS"}, \code{"Nelder-Mead"}
#'       (default: \code{"L-BFGS-B"})}
#'     \item{\code{maxit}}{Maximum number of iterations for each M-step
#'       optimization (default: 1000)}
#'     \item{\code{em_maxit}}{Maximum number of EM algorithm iterations
#'       (default: 10)}
#'     \item{\code{em_tol}}{Convergence criterion for EM algorithm based on
#'       relative change
#'       in log-likelihood (default: 1e-4)}
#'     \item{\code{tol}}{Numerical tolerance for optimization convergence
#'       (default: 1e-2)}
#'     \item{\code{verbose}}{Controls diagnostic output: \code{0}/\code{FALSE}
#'       for silent operation, \code{1}/\code{TRUE} for iteration progress,
#'       \code{2} for detailed parameter traces (default: \code{FALSE})}
#'   }
#' @param init Optional list providing initial values for model parameters.
#'   Should have the same structure as the fitted model's \code{parameters}
#'   component with elements:
#'   \describe{
#'     \item{\code{coefficients}}{A list containing:
#'       \itemize{
#'         \item \code{baseline}: Vector of B-spline coefficients for baseline
#'           hazard (length = number of spline basis functions)
#'         \item \code{hazard}: Vector of hazard parameters including
#'           association parameters (2) and survival covariates
#'         \item \code{acceleration}: Vector of longitudinal fixed effects
#'           including intercept and covariates
#'         \item \code{measurement_error_sd}: Residual standard deviation
#'           (positive scalar)
#'         \item \code{random_effect_sd}: Random effect standard deviation
#'           (positive scalar)
#'       }}
#'     \item{\code{configurations}}{Optional; if not provided, will use
#'       spline configuration from \code{spline_baseline}}
#'   }
#'   If \code{NULL}, default initial values are used (default: \code{NULL}).
#' @param parallel Logical flag enabling parallel computation for
#'   computationally intensive operations including posterior calculations,
#'   gradient evaluations, and likelihood computations. Requires \pkg{future}
#'   and \pkg{future.apply} packages
#'   (default: \code{FALSE}).
#' @param n_cores Integer specifying the number of CPU cores for parallel
#'   processing. If \code{NULL}, automatically detects and uses all available
#'   cores minus one
#'   (default: \code{NULL}).
#' @param ... Additional arguments passed to internal optimization routines.
#'
#' @return An S3 object of class \code{"JointODE"} containing fitted model
#'   results:
#'   \describe{
#'     \item{\code{parameters}}{A list containing all estimated parameters:
#'       \itemize{
#'         \item \code{coefficients}: Named list with \code{baseline} (B-spline
#'           coefficients for baseline hazard), \code{hazard} (association and
#'           survival covariate effects), \code{acceleration} (longitudinal
#'           fixed effects), \code{measurement_error_sd} (residual standard
#'           deviation), and \code{random_effect_sd} (random effect standard
#'           deviation)
#'         \item \code{configurations}: Model configuration including spline
#'           basis specifications
#'       }}
#'     \item{\code{logLik}}{Maximum log-likelihood value achieved at
#'       convergence}
#'     \item{\code{AIC}}{Akaike Information Criterion for model comparison}
#'     \item{\code{BIC}}{Bayesian Information Criterion adjusted for sample
#'       size}
#'     \item{\code{convergence}}{List containing convergence diagnostics:
#'       \itemize{
#'         \item \code{converged}: Logical indicating convergence status
#'         \item \code{em_iterations}: Number of EM iterations performed
#'         \item \code{message}: Descriptive convergence message
#'       }}
#'     \item{\code{random_effects}}{List containing random effects estimates:
#'       \itemize{
#'         \item \code{estimates}: Posterior means of subject-specific random
#'           effects
#'         \item \code{variances}: Posterior variances of random effects
#'       }}
#'     \item{\code{data}}{Processed data used for model fitting in internal
#'       format}
#'     \item{\code{control}}{List of control parameters used in optimization}
#'     \item{\code{call}}{The matched function call for reproducibility}
#'   }
#'
#' @details
#' The joint modeling framework integrates longitudinal and survival processes
#' through a shared random effects structure. The longitudinal biomarker
#' evolution is characterized by a system of ODEs that can accommodate
#' non-linear dynamics, feedback mechanisms, and complex temporal patterns.
#' The survival component employs a proportional hazards model where the
#' instantaneous risk depends on
#' features derived from the longitudinal trajectory.
#'
#' Three association structures are supported:
#' \itemize{
#'   \item Current value: hazard depends on the biomarker level at time t
#'   \item Rate of change: hazard depends on the biomarker's instantaneous
#'     slope
#'   \item Cumulative burden: hazard depends on the area under the trajectory
#'     curve
#' }
#'
#' Parameter estimation employs an Expectation-Maximization (EM) algorithm
#' with:
#' \itemize{
#'   \item E-step: Adaptive Gauss-Hermite quadrature for numerical integration
#'   \item M-step: Quasi-Newton optimization for parameter updates
#' }
#'
#' @note
#' \itemize{
#'   \item Data validation is performed automatically with informative error
#'     messages
#'   \item For high-dimensional problems, parallel computation is strongly
#'     recommended
#'   \item Convergence issues may arise with sparse event data or limited
#'     follow-up
#'   \item Initial values are computed using separate model fits when not
#'     provided
#' }
#'
#' @importFrom stats optim
#' @importFrom utils modifyList
#' @importFrom survival Surv
#' @importFrom cli cli_h2 cli_text cli_alert_success
#' @importFrom cli cli_alert_warning cli_alert_info
#' @importFrom future.apply future_lapply
#' @importFrom numDeriv jacobian
#'
#' @examples
#' \dontrun{
#' fit <- JointODE(
#'   longitudinal_formula = sim$formulas$longitudinal,
#'   longitudinal_data = sim$data$longitudinal_data,
#'   survival_formula = sim$formulas$survival,
#'   survival_data = sim$data$survival_data
#' )
#' summary(fit)
#' }
#'
#' @concept model-fitting
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
  init = NULL,
  control = list(),
  parallel = FALSE,
  n_cores = NULL,
  ...
) {
  cl <- match.call()

  # Validate and process data
  .validate(
    longitudinal_formula = longitudinal_formula,
    longitudinal_data = longitudinal_data,
    survival_formula = survival_formula,
    survival_data = survival_data,
    id = id,
    time = time,
    spline_baseline = spline_baseline,
    init = init
  )

  data_process <- .process(
    longitudinal_formula = longitudinal_formula,
    longitudinal_data = longitudinal_data,
    survival_formula = survival_formula,
    survival_data = survival_data,
    id = id,
    time = time
  )

  # Extract data dimensions
  event_times <- vapply(data_process, `[[`, numeric(1), "time")
  n_subjects <- attr(data_process, "n_subjects")
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

  # Configure splines
  spline_baseline_config <- .get_spline_config(
    x = event_times,
    degree = spline_baseline$degree,
    n_knots = spline_baseline$n_knots,
    knot_placement = spline_baseline$knot_placement,
    boundary_knots = spline_baseline$boundary_knots
  )
  spline_baseline_config$boundary_knots[1] <- 0

  # Initialize parameters
  parameters <- list(
    coefficients = list(
      baseline = rep(0, spline_baseline_config$df),
      hazard = rep(0, n_survival_covariates + 2),
      acceleration = rep(0, n_longitudinal_covariates + 3),
      measurement_error_sd = 1,
      random_effect_sd = 1
    ),
    configurations = list(baseline = spline_baseline_config)
  )

  # Override with user-provided initial values if available
  # (validation already done in .validate)
  if (!is.null(init)) {
    if (!is.null(init$coefficients)) {
      # Note: modifyList performs shallow merge, which is sufficient for
      # current parameter structure. If deeper nesting is added in future,
      # consider implementing deep merge
      parameters$coefficients <- modifyList(
        parameters$coefficients,
        init$coefficients
      )
    }

    if (!is.null(init$configurations)) {
      parameters$configurations <- modifyList(
        parameters$configurations,
        init$configurations
      )
    }
  }

  # Control settings
  control_settings <- modifyList(
    list(
      method = "L-BFGS-B",
      em_maxit = 100,
      maxit = 50,
      tol = 1e-4,
      verbose = FALSE,
      factr = 1e8
    ),
    control
  )
  em_maxit <- control_settings$em_maxit
  em_tol <- control_settings$em_tol %||% 1e-4

  # EM Algorithm setup
  converged <- FALSE
  old_loglik <- -Inf
  n_baseline <- spline_baseline_config$df
  n_hazard <- n_survival_covariates + 2

  # Convert verbose to numeric level
  verbose_level <- if (is.logical(control_settings$verbose)) {
    as.numeric(control_settings$verbose)
  } else {
    control_settings$verbose
  }

  if (verbose_level > 0) {
    cli::cli_h2("EM Algorithm Optimization")
    cli::cli_alert_info(
      "Settings: {em_maxit} iterations | tolerance: {em_tol}"
    )
    cli::cli_text("")
  }

  for (em_iter in seq_len(em_maxit)) {
    # E-step
    posteriors <- .compute_posteriors(
      data_process,
      parameters,
      parallel = parallel,
      n_cores = n_cores,
      return_ode_solutions = TRUE
    )

    # Update variance components
    variance_update <- .update_variance_components(
      data_process,
      parameters,
      posteriors,
      posteriors$ode_solutions
    )
    parameters$coefficients$measurement_error_sd <-
      variance_update$measurement_error_sd
    parameters$coefficients$random_effect_sd <-
      variance_update$random_effect_sd

    fixed_parameters <- list(
      measurement_error_sd = parameters$coefficients$measurement_error_sd,
      random_effect_sd = parameters$coefficients$random_effect_sd
    )

    # M-step
    res <- optim(
      par = unlist(
        parameters$coefficients[c(
          "baseline",
          "hazard",
          "acceleration"
        )],
        use.names = FALSE
      ),
      fn = .compute_objective_joint,
      gr = .compute_gradient_joint,
      method = control_settings$method,
      control = list(
        pgtol = control_settings$tol,
        maxit = control_settings$maxit,
        factr = control_settings$factr,
        trace = if (verbose_level >= 3) 1 else 0
      ),
      data_list = data_process,
      posteriors = posteriors,
      configurations = list(baseline = spline_baseline_config),
      fixed_parameters = fixed_parameters,
      parallel = parallel,
      n_cores = n_cores
    )

    # Update parameters
    idx_end <- cumsum(c(
      n_baseline,
      n_hazard,
      length(res$par) - n_baseline - n_hazard
    ))
    idx_start <- c(1, idx_end[-length(idx_end)] + 1)
    parameters$coefficients$baseline <- res$par[idx_start[1]:idx_end[1]]
    parameters$coefficients$hazard <- res$par[idx_start[2]:idx_end[2]]
    parameters$coefficients$acceleration <- res$par[idx_start[3]:idx_end[3]]

    # Track convergence
    new_loglik <- -res$value
    loglik_change <- new_loglik - old_loglik

    if (verbose_level > 0) {
      cli::cli_text(
        "[{sprintf('%3d', em_iter)}/{em_maxit}] ",
        "LogLik: {sprintf('%.4f', new_loglik)} ",
        "(Change: {sprintf('%.4f', loglik_change)}) ",
        "\u03c3_e={sprintf('%.3f',
                     parameters$coefficients$measurement_error_sd)} ",
        "\u03c3_b={sprintf('%.3f', parameters$coefficients$random_effect_sd)}"
      )

      # Detailed parameter output for verbose=2
      if (verbose_level >= 2) {
        cli::cli_text(
          "  \u03b7: [{paste(sprintf('%.3f',
            head(parameters$coefficients$baseline, 5)), collapse=', ')}",
          if (length(parameters$coefficients$baseline) > 5) "..." else "",
          "]"
        )
        cli::cli_text(
          "  \u03b1: [{paste(sprintf('%.3f',
            parameters$coefficients$hazard[1:2]), collapse=', ')}]"
        )
        if (n_hazard > 2) {
          cli::cli_text(
            "  \u03c6: [{paste(sprintf('%.3f',
              parameters$coefficients$hazard[3:n_hazard]), collapse=', ')}]"
          )
        }
        cli::cli_text(
          "  \u03b2: [{paste(sprintf('%.3f',
            parameters$coefficients$acceleration), collapse=', ')}]"
        )
        cli::cli_text("")
      }
    }

    # Check convergence
    if (abs(loglik_change) < em_tol && em_iter > 1) {
      converged <- TRUE
      if (verbose_level > 0) {
        cli::cli_text("")
        cli::cli_alert_success(
          "Converged after {em_iter} iterations ",
          "(LogLik: {sprintf('%.4f', new_loglik)})"
        )
      }
      break
    }

    old_loglik <- new_loglik
  }

  if (!converged && verbose_level > 0) {
    cli::cli_text("")
    cli::cli_alert_warning(
      "Did not converge within {em_maxit} iterations ",
      "(LogLik: {sprintf('%.4f', new_loglik)})"
    )
    cli::cli_alert_info(
      "Try increasing em_maxit or adjusting initial values"
    )
  }

  # Final computations
  final_posteriors <- .compute_posteriors(
    data_process,
    parameters,
    parallel = parallel,
    n_cores = n_cores,
    return_ode_solutions = TRUE
  )

  log_lik <- -res$value
  n_params <- length(res$par) + 2
  aic <- -2 * log_lik + 2 * n_params
  bic <- -2 * log_lik + n_params * log(n_subjects)

  # Compute variance-covariance matrix
  vcov_matrix <- if (!converged) {
    if (verbose_level > 0) {
      cli::cli_alert_warning(
        "Skipping variance-covariance computation (EM did not converge)"
      )
    }
    NULL
  } else {
    tryCatch(
      {
        if (verbose_level > 0) {
          cli::cli_alert_info("Computing variance-covariance matrix...")
        }

        # Compute Hessian via numerical differentiation of gradient
        # H = d²(-logL)/dθ² ≈ J(∇(-logL))
        H <- numDeriv::jacobian(
          func = function(p) {
            .compute_gradient_joint(
              p,
              data_list = data_process,
              posteriors = final_posteriors,
              configurations = list(baseline = spline_baseline_config),
              fixed_parameters = fixed_parameters,
              parallel = parallel,
              n_cores = n_cores
            )
          },
          x = res$par,
          method = "Richardson" # More accurate for second derivatives
        )

        # Ensure symmetry (numerical errors can break symmetry)
        H <- 0.5 * (H + t(H))

        # Check condition before inversion
        eigen_vals <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
        if (min(eigen_vals) <= 0) {
          stop("Hessian is not positive definite")
        }

        # Invert to get covariance matrix
        V <- solve(H)

        # Set parameter names
        n_accel <- length(res$par) - n_baseline - n_hazard
        param_names <- c(
          paste0("baseline:", seq_len(n_baseline)),
          paste0("hazard:", c("alpha1", "alpha2")),
          if (n_hazard > 2) paste0("hazard:phi", seq_len(n_hazard - 2)),
          paste0("longitudinal:beta", seq_len(n_accel))
        )
        dimnames(V) <- list(param_names, param_names)

        if (verbose_level > 0) {
          cli::cli_alert_success("Variance-covariance matrix computed")
        }
        V
      },
      error = function(e) {
        if (verbose_level > 0) {
          cli::cli_alert_warning(
            "Variance-covariance matrix computation failed: {e$message}"
          )
          cli::cli_alert_info(
            "Standard errors will not be available in summary output"
          )
        }
        NULL
      }
    )
  }

  # Return fitted model
  structure(
    list(
      parameters = parameters,
      logLik = log_lik,
      AIC = aic,
      BIC = bic,
      convergence = list(
        converged = converged,
        em_iterations = if (converged) em_iter else em_maxit,
        message = sprintf(
          "EM algorithm %s after %d iterations",
          if (converged) "converged" else "did not converge within",
          if (converged) em_iter else em_maxit
        )
      ),
      random_effects = list(
        estimates = final_posteriors$b,
        variances = final_posteriors$v
      ),
      vcov = vcov_matrix,
      data = data_process,
      control = control_settings,
      call = cl
    ),
    class = "JointODE"
  )
}

#' Summary of JointODE Model
#'
#' @param object A JointODE object
#' @param ... Additional arguments
#' @return A summary.JointODE object with coefficients and test statistics
#'
#' @concept model-summary
#' @importFrom stats coef pnorm
#' @export
summary.JointODE <- function(object, ...) {
  coefs <- coef(object)
  se <- if (!is.null(object$vcov)) {
    sqrt(diag(object$vcov))
  } else {
    rep(NA_real_, length(coefs))
  }

  structure(
    list(
      call = object$call,
      coefficients = cbind(
        Estimate = coefs,
        `Std. Error` = se,
        `z value` = coefs / se,
        `Pr(>|z|)` = 2 * pnorm(-abs(coefs / se))
      ),
      sigma = with(
        object$parameters$coefficients,
        c(sigma_e = measurement_error_sd, sigma_b = random_effect_sd)
      ),
      logLik = object$logLik,
      AIC = object$AIC,
      BIC = object$BIC,
      nobs = attr(object$data, "n_subjects"),
      convergence = object$convergence
    ),
    class = "summary.JointODE"
  )
}

#' @concept model-summary
#' @importFrom stats printCoefmat
#' @export
print.summary.JointODE <- function(
  x,
  digits = max(3L, getOption("digits") - 3L),
  signif.stars = getOption("show.signif.stars"),
  ...
) {
  cat("\nCall:\n")
  print(x$call)
  cat("\nVariance components:\n")
  print(x$sigma, digits = digits)
  cat("\nFixed effects:\n")
  printCoefmat(
    x$coefficients,
    digits = digits,
    signif.stars = signif.stars,
    ...
  )
  cat("\n---")
  cat("\nLog-likelihood:", x$logLik, "  AIC:", x$AIC, "  BIC:", x$BIC)
  cat("\nN =", x$nobs, " Convergence:", x$convergence$message, "\n")
  invisible(x)
}

#' Extract Model Coefficients
#'
#' @param object A JointODE object
#' @param ... Additional arguments
#' @return Named numeric vector of fixed effects coefficients
#' @concept model-inspection
#' @export
coef.JointODE <- function(object, ...) {
  cf <- object$parameters$coefficients
  coefs <- c(cf$baseline, cf$hazard, cf$acceleration)
  names(coefs) <- c(
    paste0("baseline:", seq_along(cf$baseline)),
    paste0(
      "hazard:",
      c(
        "alpha1",
        "alpha2",
        if (length(cf$hazard) > 2) paste0("phi", seq_len(length(cf$hazard) - 2))
      )
    ),
    paste0("longitudinal:beta", seq_along(cf$acceleration))
  )
  coefs
}

#' Extract Variance-Covariance Matrix
#'
#' @param object A JointODE object
#' @param ... Additional arguments
#' @return Variance-covariance matrix of fixed effects
#' @concept model-inspection
#' @export
vcov.JointODE <- function(object, ...) {
  object$vcov
}

#' Extract Log-Likelihood
#'
#' @param object A JointODE object
#' @param ... Additional arguments
#' @return Log-likelihood with df and nobs attributes
#' @concept model-inspection
#' @importFrom stats coef
#' @export
logLik.JointODE <- function(object, ...) {
  structure(
    object$logLik,
    df = length(coef(object)) + 2,
    nobs = attr(object$data, "n_subjects"),
    class = "logLik"
  )
}

#' Print JointODE Model
#'
#' @param x A JointODE object
#' @param digits Number of digits for numeric output
#' @param ... Additional arguments
#' @return Invisibly returns the object
#' @concept model-display
#' @importFrom stats coef
#' @export
print.JointODE <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nJoint ODE Model\n")
  cat("Call: ")
  print(x$call)
  cat(
    "\nLog-likelihood:",
    format(x$logLik, digits = digits),
    "on",
    length(coef(x)) + 2,
    "degrees of freedom\n"
  )
  cat(
    "AIC:",
    format(x$AIC, digits = digits),
    "  BIC:",
    format(x$BIC, digits = digits),
    "\n"
  )
  invisible(x)
}

#' Predict Method for JointODE Models
#'
#' @description
#' Computes predictions for biomarker trajectories, velocities, accelerations,
#' and survival functions from a fitted JointODE model. This method uses the
#' fitted model parameters and random effects to solve the ODE system and
#' generate predictions for each subject.
#'
#' @param object A fitted JointODE model object.
#' @param times Optional numeric vector of time points at which to evaluate
#'   predictions. If NULL, uses the union of observation times and event time
#'   for each subject.
#' @param parallel Logical flag for parallel computation. Default is FALSE.
#' @param n_cores Number of CPU cores for parallel processing. If NULL,
#'   uses all available cores minus one.
#' @param ... Additional arguments (currently unused).
#'
#' @return A list containing predictions for each subject, where each
#'   element contains:
#'   \describe{
#'     \item{\code{id}}{Subject identifier}
#'     \item{\code{times}}{Time points for predictions}
#'     \item{\code{biomarker}}{Predicted biomarker values with random effects}
#'     \item{\code{velocity}}{Predicted biomarker velocity}
#'     \item{\code{acceleration}}{Predicted biomarker acceleration}
#'     \item{\code{cumhazard}}{Cumulative hazard values}
#'     \item{\code{survival}}{Survival probability values}
#'   }
#'
#' @details
#' The prediction process involves:
#' \enumerate{
#'   \item Using the fitted model parameters and random effects
#'   \item Solving the ODE system for each subject at specified time points
#'   \item Computing biomarker trajectories with random effects included
#'   \item Calculating cumulative hazard and survival probabilities
#' }
#'
#' Random effects from the fitted model are always included in predictions
#' to provide subject-specific trajectories.
#' @concept model-prediction
#' @importFrom stats predict
#' @export
predict.JointODE <- function(
  object,
  times = NULL,
  parallel = FALSE,
  n_cores = NULL,
  ...
) {
  # Extract model components
  parameters <- object$parameters
  data_process <- object$data
  random_effects <- object$random_effects$estimates

  # Function to compute predictions for a single subject
  compute_subject_predictions <- function(i) {
    subject_data <- data_process[[i]]
    b_i <- random_effects[i]

    # Determine time points for prediction
    if (!is.null(times)) {
      pred_times <- sort(unique(times))
    } else {
      # Use observation times (without event time for cleaner output)
      pred_times <- sort(unique(
        subject_data$longitudinal$times,
        subject_data$time
      ))
    }

    # Solve ODE with specified times
    ode_solution <- .solve_joint_ode(
      subject_data,
      parameters,
      sensitivity_type = "basic",
      times = pred_times
    )

    cumhazard_values <- ode_solution$cum_hazard * exp(b_i)
    survival_values <- exp(-cumhazard_values)

    predictions <- list(
      id = subject_data$id,
      times = pred_times,
      biomarker = ode_solution$biomarker + b_i,
      velocity = ode_solution$velocity,
      acceleration = ode_solution$acceleration,
      survival = survival_values
    )

    predictions
  }

  # Compute predictions for all subjects
  n_subjects <- length(data_process)

  if (parallel) {
    cleanup <- .setup_parallel_plan(n_cores)
    on.exit(cleanup(), add = TRUE)

    prediction_results <- future.apply::future_lapply(
      seq_len(n_subjects),
      compute_subject_predictions,
      future.seed = TRUE
    )
  } else {
    prediction_results <- lapply(
      seq_len(n_subjects),
      compute_subject_predictions
    )
  }

  # Return predictions
  prediction_results
}
