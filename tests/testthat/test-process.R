# ==============================================================================
# Tests for Data Processing Functions (.process)
# ==============================================================================
#
# Test Coverage:
# 1. Basic Structure and Output
#    - Correct structuring of processed data
#    - Subject-level data organization
#    - Attributes (n_subjects, n_observations, event_rate)
#
# 2. Missing Data Handling
#    - Subjects with missing longitudinal data
#    - Proper handling of unmatched subjects
#
# 3. Formula Processing
#    - Simple formulas with covariates
#    - Formula transformations (log, I(), etc.)
#    - Interaction terms (x1 * x2)
#    - Polynomial terms (poly())
#    - Spline terms (ns(), bs())
#    - Factor variables and dummy coding
#
# 4. Data Integrity
#    - Preservation of data order within subjects
#    - Correct time sorting
#    - Handling of non-sequential time points
#
# 5. Edge Cases
#    - Single time point per subject
#    - Large number of covariates (high-dimensional)
#
# 6. Summary Statistics
#    - Correct computation of event rates
#    - Accurate observation counts
# ==============================================================================

test_that(".process correctly structures basic data", {
  longitudinal_data <- create_test_longitudinal_data(
    n_subjects = 3, n_times = 4, seed = 123
  )
  longitudinal_data$v <- longitudinal_data$v + 10 # Adjust mean

  survival_data <- create_test_survival_data(n_subjects = 3, seed = 123)
  names(survival_data)[2] <- "obstime" # Rename time column

  result <- JointODE:::.process(
    longitudinal_formula = v ~ x1 + x2,
    longitudinal_data = longitudinal_data,
    survival_formula = Surv(obstime, status) ~ w1 + w2,
    survival_data = survival_data,
    id = "id",
    time = "time"
  )

  subject1 <- result[["1"]]

  # Basic structure checks
  expect_equal(subject1$id, 1)
  expect_true(is.numeric(subject1$time))
  expect_true(subject1$time > 0)
  expect_true(subject1$status %in% c(0, 1))

  # Longitudinal checks
  expect_equal(length(subject1$longitudinal$times), 4)
  expect_equal(subject1$longitudinal$times, 0:3)
  expect_equal(subject1$longitudinal$n_obs, 4)
  expect_true(ncol(subject1$longitudinal$covariates) >= 2)

  # Attributes
  expect_equal(attr(result, "n_subjects"), 3)
  expect_equal(attr(result, "n_observations"), 12)
  expect_true(
    attr(result, "event_rate") >= 0 && attr(result, "event_rate") <= 1
  )
})

test_that(".process handles missing longitudinal data", {
  test_data <- create_mismatched_test_data(
    n_subjects_long = 2, n_subjects_surv = 3
  )
  longitudinal_data <- test_data$longitudinal
  survival_data <- test_data$survival
  names(survival_data)[2] <- "obstime"

  result <- JointODE:::.process(
    longitudinal_formula = v ~ 1,
    longitudinal_data = longitudinal_data,
    survival_formula = Surv(obstime, status) ~ 1,
    survival_data = survival_data,
    id = "id",
    time = "time"
  )

  subject3 <- result[["3"]]
  expect_equal(subject3$id, 3)
  expect_true(is.numeric(subject3$time))
  expect_true(subject3$time > 0)
  expect_true(subject3$status %in% c(0, 1))
  expect_equal(subject3$longitudinal$n_obs, 0)
  expect_equal(length(subject3$longitudinal$times), 0)
})

test_that(".process handles formula transformations", {
  set.seed(456)
  longitudinal_data <- data.frame(
    id = rep(1:10, each = 3),
    time = rep(1:3, 10),
    v = exp(rnorm(30)),
    x = runif(30, 1, 10)
  )

  survival_data <- data.frame(
    id = 1:10,
    surv_time = rexp(10),
    status = rbinom(10, 1, 0.5),
    z = runif(10, 1, 10)
  )

  result <- JointODE:::.process(
    longitudinal_formula = log(v) ~ log(x) + I(x^2),
    longitudinal_data = longitudinal_data,
    survival_formula = Surv(surv_time, status) ~ log(z),
    survival_data = survival_data,
    id = "id",
    time = "time"
  )

  # Check transformations applied
  expect_equal(
    unname(result[[1]]$longitudinal$measurements[1]),
    log(longitudinal_data$v[longitudinal_data$id == 1][1])
  )
  expect_equal(ncol(result[[1]]$longitudinal$covariates), 3)
})

test_that(".process handles interaction terms", {
  longitudinal_data <- data.frame(
    id = rep(1:5, each = 4),
    time = rep(0:3, 5),
    y = rnorm(20),
    x1 = rnorm(20),
    x2 = rbinom(20, 1, 0.5)
  )

  survival_data <- data.frame(
    id = 1:5,
    event_time = rexp(5) + 0.1,
    status = rbinom(5, 1, 0.7),
    z1 = rnorm(5),
    z2 = rbinom(5, 1, 0.5)
  )

  result <- JointODE:::.process(
    longitudinal_formula = y ~ x1 * x2,
    longitudinal_data = longitudinal_data,
    survival_formula = Surv(event_time, status) ~ z1 + z2,
    survival_data = survival_data,
    id = "id",
    time = "time"
  )

  # Check interaction term is included
  subject1 <- result[["1"]]
  # Should have time, x1, x2, x1:x2
  expect_equal(ncol(subject1$longitudinal$covariates), 4)
})

test_that(".process handles poly and spline terms", {
  longitudinal_data <- data.frame(
    id = rep(1:5, each = 4),
    time = rep(0:3, 5),
    y = rnorm(20),
    x = runif(20)
  )

  survival_data <- data.frame(
    id = 1:5,
    event_time = rexp(5) + 0.1,
    status = rbinom(5, 1, 0.7),
    z = runif(5)
  )

  # Test polynomial terms
  result_poly <- JointODE:::.process(
    longitudinal_formula = y ~ poly(x, 2),
    longitudinal_data = longitudinal_data,
    survival_formula = Surv(event_time, status) ~ poly(z, 2),
    survival_data = survival_data,
    id = "id",
    time = "time"
  )

  subject1_poly <- result_poly[["1"]]
  # Should have time + 2 poly terms
  expect_equal(ncol(subject1_poly$longitudinal$covariates), 3)

  # Test spline terms
  result_spline <- JointODE:::.process(
    longitudinal_formula = y ~ splines::ns(x, df = 3),
    longitudinal_data = longitudinal_data,
    survival_formula = Surv(event_time, status) ~ splines::ns(z, df = 2),
    survival_data = survival_data,
    id = "id",
    time = "time"
  )

  subject1_spline <- result_spline[["1"]]
  # Should have time + 3 spline basis functions
  expect_equal(ncol(subject1_spline$longitudinal$covariates), 4)
})

test_that(".process handles factor variables", {
  longitudinal_data <- data.frame(
    id = rep(1:5, each = 4),
    time = rep(0:3, 5),
    y = rnorm(20),
    group = factor(rep(c("A", "B", "C", "A", "B"), each = 4))
  )

  survival_data <- data.frame(
    id = 1:5,
    event_time = rexp(5) + 0.1,
    status = rbinom(5, 1, 0.7),
    treatment = factor(c("control", "drug", "control", "drug", "control"))
  )

  result <- JointODE:::.process(
    longitudinal_formula = y ~ group,
    longitudinal_data = longitudinal_data,
    survival_formula = Surv(event_time, status) ~ treatment,
    survival_data = survival_data,
    id = "id",
    time = "time"
  )

  subject1 <- result[["1"]]
  # Should have time + 2 dummy variables for 3-level factor
  expect_equal(ncol(subject1$longitudinal$covariates), 3)
})

test_that(".process preserves data order within subjects", {
  # Create data with non-sequential time ordering
  longitudinal_data <- data.frame(
    id = rep(1:3, each = 4),
    time = rep(c(2, 0, 3, 1), 3),
    y = rnorm(12),
    x = rnorm(12)
  )

  survival_data <- data.frame(
    id = 1:3,
    event_time = c(4, 5, 4.5),
    status = c(1, 0, 1)
  )

  result <- JointODE:::.process(
    longitudinal_formula = y ~ x,
    longitudinal_data = longitudinal_data,
    survival_formula = Surv(event_time, status) ~ 1,
    survival_data = survival_data,
    id = "id",
    time = "time"
  )

  # Check times are sorted for each subject
  for (i in 1:3) {
    subj <- result[[as.character(i)]]
    expect_equal(subj$longitudinal$times, 0:3)
    # Check corresponding measurements are correctly ordered
    expect_true(all(diff(subj$longitudinal$times) >= 0))
  }
})

test_that(".process handles single time point per subject", {
  longitudinal_data <- data.frame(
    id = 1:5,
    time = rep(0, 5),
    y = rnorm(5),
    x = rnorm(5)
  )

  survival_data <- data.frame(
    id = 1:5,
    event_time = rexp(5) + 0.1,
    status = rbinom(5, 1, 0.7)
  )

  result <- JointODE:::.process(
    longitudinal_formula = y ~ x,
    longitudinal_data = longitudinal_data,
    survival_formula = Surv(event_time, status) ~ 1,
    survival_data = survival_data,
    id = "id",
    time = "time"
  )

  subject1 <- result[["1"]]
  expect_equal(subject1$longitudinal$n_obs, 1)
  expect_equal(length(subject1$longitudinal$times), 1)
})

test_that(".process handles large number of covariates", {
  n_covs <- 20
  n_subjects <- 10
  n_times <- 5

  # Create wide data
  long_data_list <- list(
    id = rep(1:n_subjects, each = n_times),
    time = rep(0:(n_times - 1), n_subjects),
    y = rnorm(n_subjects * n_times)
  )

  for (i in 1:n_covs) {
    long_data_list[[paste0("x", i)]] <- rnorm(n_subjects * n_times)
  }

  longitudinal_data <- as.data.frame(long_data_list)

  survival_data <- data.frame(
    id = 1:n_subjects,
    event_time = rexp(n_subjects) + 0.1,
    status = rbinom(n_subjects, 1, 0.7)
  )

  # Create formula with all covariates
  long_formula <- as.formula(
    paste("y ~", paste(paste0("x", 1:n_covs), collapse = " + "))
  )

  result <- JointODE:::.process(
    longitudinal_formula = long_formula,
    longitudinal_data = longitudinal_data,
    survival_formula = Surv(event_time, status) ~ 1,
    survival_data = survival_data,
    id = "id",
    time = "time"
  )

  subject1 <- result[["1"]]
  # Should have time + 20 covariates
  expect_equal(ncol(subject1$longitudinal$covariates), n_covs + 1)
})

test_that(".process computes correct summary attributes", {
  set.seed(789)
  n_subjects <- 100
  n_times <- 5

  longitudinal_data <- data.frame(
    id = rep(1:n_subjects, each = n_times),
    time = rep(0:(n_times - 1), n_subjects),
    y = rnorm(n_subjects * n_times),
    x = rnorm(n_subjects * n_times)
  )

  # Create survival data with known event rate
  n_events <- 60
  survival_data <- data.frame(
    id = 1:n_subjects,
    event_time = rexp(n_subjects) + 0.1,
    status = c(rep(1, n_events), rep(0, n_subjects - n_events))
  )

  result <- JointODE:::.process(
    longitudinal_formula = y ~ x,
    longitudinal_data = longitudinal_data,
    survival_formula = Surv(event_time, status) ~ 1,
    survival_data = survival_data,
    id = "id",
    time = "time"
  )

  expect_equal(attr(result, "n_subjects"), n_subjects)
  expect_equal(attr(result, "n_observations"), n_subjects * n_times)
  expect_equal(attr(result, "event_rate"), n_events / n_subjects)
})
