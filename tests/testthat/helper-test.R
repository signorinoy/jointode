# Helper functions for creating test data

#' Create standard longitudinal test data
#' @param n_subjects Number of subjects
#' @param n_times Number of time points per subject
#' @param seed Random seed for reproducibility
create_test_longitudinal_data <- function(n_subjects = 3, n_times = 3,
                                          seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  data.frame(
    id = rep(seq_len(n_subjects), each = n_times),
    time = rep(seq(0, n_times - 1), n_subjects),
    v = rnorm(n_subjects * n_times),
    x1 = rnorm(n_subjects * n_times),
    x2 = factor(rep(seq_len(n_subjects), each = n_times))
  )
}

#' Create standard survival test data
#' @param n_subjects Number of subjects
#' @param event_rate Proportion of events (vs censored)
#' @param seed Random seed for reproducibility
create_test_survival_data <- function(n_subjects = 3, event_rate = 0.67,
                                      seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  data.frame(
    id = seq_len(n_subjects),
    time = rexp(n_subjects, rate = 0.3) + 2,
    status = rbinom(n_subjects, 1, event_rate),
    w1 = rnorm(n_subjects),
    w2 = factor(sample(c("low", "medium", "high"), n_subjects, replace = TRUE))
  )
}

#' Create minimal test data with specified characteristics
#' @param n_subjects Number of subjects
#' @param n_times Number of time points per subject
#' @param with_covariates Whether to include covariates
create_minimal_test_data <- function(n_subjects = 2, n_times = 2,
                                     with_covariates = FALSE) {
  longitudinal_data <- data.frame(
    id = rep(seq_len(n_subjects), each = n_times),
    time = rep(seq(0, n_times - 1), n_subjects),
    v = rnorm(n_subjects * n_times)
  )

  survival_data <- data.frame(
    id = seq_len(n_subjects),
    time = runif(n_subjects, min = n_times, max = n_times + 2),
    status = sample(0:1, n_subjects, replace = TRUE)
  )

  if (with_covariates) {
    longitudinal_data$x1 <- rnorm(nrow(longitudinal_data))
    survival_data$w1 <- rnorm(nrow(survival_data))
  }

  list(longitudinal = longitudinal_data, survival = survival_data)
}

#' Create test data with missing subjects
#' @param n_subjects_long Number of subjects in longitudinal data
#' @param n_subjects_surv Number of subjects in survival data
#' @param n_times Number of time points per subject
create_mismatched_test_data <- function(n_subjects_long, n_subjects_surv,
                                        n_times = 3) {
  longitudinal_data <- if (n_subjects_long > 0) {
    data.frame(
      id = rep(seq_len(n_subjects_long), each = n_times),
      time = rep(seq(0, n_times - 1), n_subjects_long),
      v = rnorm(n_subjects_long * n_times)
    )
  } else {
    data.frame(id = numeric(0), time = numeric(0), v = numeric(0))
  }

  survival_data <- if (n_subjects_surv > 0) {
    data.frame(
      id = seq_len(n_subjects_surv),
      time = rexp(n_subjects_surv, 0.3) + 2,
      status = rbinom(n_subjects_surv, 1, 0.5)
    )
  } else {
    data.frame(id = numeric(0), time = numeric(0), status = numeric(0))
  }

  list(longitudinal = longitudinal_data, survival = survival_data)
}


estimate_bspline_coef <- function(x, f0, config) {
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

  # Fit coefficients using least squares
  spline_coefficients <- solve(
    t(basis_matrix) %*% basis_matrix,
    t(basis_matrix) %*% y_target
  )
  as.vector(spline_coefficients)
}

setup_test_environment <- function() {
  # Define default configurations
  spline_baseline <- list(
    degree = 3,
    n_knots = 5,
    knot_placement = "quantile",
    boundary_knots = NULL
  )

  spline_index <- list(
    degree = 3,
    n_knots = 10,
    knot_placement = "quantile",
    boundary_knots = NULL
  )

  # Define functions once
  g_0 <- function(u) 2 * tanh(u / 3)
  lambda_0 <- function(t) log((1.5 / 8) * (t / 8)^0.5)

  # Generate and process data
  data <- simulate()
  data_process <- .process(
    longitudinal_data = data$longitudinal_data,
    longitudinal_formula = v ~ x1 + x2,
    survival_data = data$survival_data,
    survival_formula = Surv(time, status) ~ w1 + w2,
    id = "id",
    time = "time"
  )

  # Use override or computed number of subjects
  n_subjects <- length(data_process)

  # Extract times once
  times <- data$survival_data[, "time"]

  # Compute baseline spline coefficients
  baseline_spline_coefficients <- estimate_bspline_coef(
    times, lambda_0, spline_baseline
  )

  # Define coefficients
  hazard_coefficients <- c(0.5, 0.4, -0.5, 0.2, -0.15)
  index_beta_raw <- c(-0.3, -0.5, 0, 0.2, 0.1, 0.05)
  index_coefficients <- index_beta_raw / sqrt(sum(index_beta_raw^2))

  # Compute index variable efficiently
  z_cols <- c("biomarker", "velocity", "x1", "x2", "time")
  z_coefs <- c(-0.3, -0.5, 0.2, 0.1, 0.05)
  Z <- as.matrix(data$longitudinal_data[, z_cols])
  u <- drop(Z %*% z_coefs) # drop() removes unnecessary dimensions

  # Compute index spline coefficients
  index_spline_coefficients <- estimate_bspline_coef(
    u, g_0, spline_index
  )

  # Create spline configurations
  spline_baseline_config <- .get_spline_config(
    x = times,
    degree = spline_baseline$degree,
    n_knots = spline_baseline$n_knots,
    knot_placement = spline_baseline$knot_placement,
    boundary_knots = spline_baseline$boundary_knots
  )

  spline_index_config <- .get_spline_config(
    x = u,
    degree = spline_index$degree,
    n_knots = spline_index$n_knots,
    knot_placement = spline_index$knot_placement,
    boundary_knots = spline_index$boundary_knots
  )

  # Generate posteriors
  posteriors <- list(
    b = rnorm(n_subjects, 0, 0.2),
    v = runif(n_subjects, 0.01, 0.1),
    exp_b = exp(rnorm(n_subjects, 0, 0.2))
  )

  # Return organized structure
  list(
    data_list = data_process,
    config = list(
      baseline = spline_baseline_config,
      index = spline_index_config
    ),
    params = list(
      baseline = baseline_spline_coefficients,
      hazard = hazard_coefficients,
      index_g = index_spline_coefficients,
      index_beta = index_coefficients,
      measurement_error_sd = 0.1,  # Changed from 1e-2 to avoid numerical issues
      random_effect_sd = 0.1       # Changed from 1e-2 to avoid numerical issues
    ),
    posteriors = posteriors
  )
}
