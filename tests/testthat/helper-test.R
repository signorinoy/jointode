# Helper functions for creating test data

#' Create standard longitudinal test data
#' @param n_subjects Number of subjects
#' @param n_times Number of time points per subject
#' @param seed Random seed for reproducibility
create_test_longitudinal_data <- function(
  n_subjects = 3,
  n_times = 3,
  seed = NULL
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

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
create_test_survival_data <- function(
  n_subjects = 3,
  event_rate = 0.67,
  seed = NULL
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

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
create_minimal_test_data <- function(
  n_subjects = 2,
  n_times = 2,
  with_covariates = FALSE
) {
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
create_mismatched_test_data <- function(
  n_subjects_long,
  n_subjects_surv,
  n_times = 3
) {
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

# ==============================================================================
# Helper Functions for Gradient Tests
# ==============================================================================

#' Compare gradient components with detailed diagnostics
#' @param analytical Analytical gradient vector
#' @param numerical Numerical gradient vector
#' @param name Component name for error messages
#' @param tolerance Comparison tolerance (default: 1e-4)
compare_gradient_component <- function(
  analytical,
  numerical,
  name,
  tolerance = 1e-2
) {
  rel_error <- max(abs((analytical - numerical) / (abs(numerical) + 1e-10)))

  expect_lt(rel_error, tolerance, label = paste0(name, " gradient"))
}

#' Load and process test data from sim dataset
#' @param n_subjects Number of subjects to use (NULL = all)
#' @return List with data, parameters, and n_subjects
load_test_data <- function(n_subjects = NULL) {
  # Process data using sim dataset
  data_processed <- .process(
    longitudinal_data = sim$data$longitudinal_data,
    longitudinal_formula = sim$formulas$longitudinal,
    survival_data = sim$data$survival_data,
    survival_formula = sim$formulas$survival,
    id = "id",
    time = "time"
  )

  parameters <- sim$parameters

  # Subset if requested for faster testing
  if (!is.null(n_subjects)) {
    n_subjects <- min(n_subjects, length(data_processed))
    data_processed <- data_processed[1:n_subjects]
  }

  list(
    data = data_processed,
    parameters = parameters
  )
}
