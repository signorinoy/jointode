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
  data.long <- data.frame(
    id = rep(seq_len(n_subjects), each = n_times),
    time = rep(seq(0, n_times - 1), n_subjects),
    v = rnorm(n_subjects * n_times)
  )

  data.surv <- data.frame(
    id = seq_len(n_subjects),
    time = runif(n_subjects, min = n_times, max = n_times + 2),
    status = sample(0:1, n_subjects, replace = TRUE)
  )

  if (with_covariates) {
    data.long$x1 <- rnorm(nrow(data.long))
    data.surv$w1 <- rnorm(nrow(data.surv))
  }

  list(longitudinal = data.long, survival = data.surv)
}

#' Create test data with missing subjects
#' @param n_subjects_long Number of subjects in longitudinal data
#' @param n_subjects_surv Number of subjects in survival data
#' @param n_times Number of time points per subject
create_mismatched_test_data <- function(n_subjects_long, n_subjects_surv,
                                        n_times = 3) {
  data.long <- if (n_subjects_long > 0) {
    data.frame(
      id = rep(seq_len(n_subjects_long), each = n_times),
      time = rep(seq(0, n_times - 1), n_subjects_long),
      v = rnorm(n_subjects_long * n_times)
    )
  } else {
    data.frame(id = numeric(0), time = numeric(0), v = numeric(0))
  }

  data.surv <- if (n_subjects_surv > 0) {
    data.frame(
      id = seq_len(n_subjects_surv),
      time = rexp(n_subjects_surv, 0.3) + 2,
      status = rbinom(n_subjects_surv, 1, 0.5)
    )
  } else {
    data.frame(id = numeric(0), time = numeric(0), status = numeric(0))
  }

  list(longitudinal = data.long, survival = data.surv)
}
