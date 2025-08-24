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
