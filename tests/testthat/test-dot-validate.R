test_that(".validate() checks input types correctly", {
  # Create valid test data
  longitudinal_data <- create_test_longitudinal_data(
    n_subjects = 3, n_times = 3, seed = 123
  )[, c("id", "time", "v", "x1")]
  survival_data <- create_test_survival_data(
    n_subjects = 3, seed = 123
  )[, c("id", "time", "status")]
  survival_data$w1 <- rnorm(3)

  longitudinal_formula <- v ~ x1
  survival_formula <- Surv(time, status) ~ w1

  # Test invalid formula types
  expect_error(
    .validate(
      "not a formula", longitudinal_data, survival_formula, survival_data,
      "id", "time"
    ),
    "longitudinal_formula must be a formula"
  )

  expect_error(
    .validate(
      longitudinal_formula, longitudinal_data, "not a formula", survival_data,
      "id", "time"
    ),
    "survival_formula must be a formula"
  )

  # Test invalid data frame types
  expect_error(
    .validate(
      longitudinal_formula, "not a df", survival_formula, survival_data,
      "id", "time"
    ),
    "longitudinal_data must be a data frame"
  )

  expect_error(
    .validate(
      longitudinal_formula, longitudinal_data, survival_formula, list(),
      "id", "time"
    ),
    "survival_data must be a data frame"
  )

  # Test invalid ID and time types
  expect_error(
    .validate(
      longitudinal_formula, longitudinal_data, survival_formula, survival_data,
      123, "time"
    ),
    "id must be a single character string"
  )

  expect_error(
    .validate(
      longitudinal_formula, longitudinal_data, survival_formula, survival_data,
      "id", c("t1", "t2")
    ),
    "time must be a single character string"
  )
})

test_that(".validate() checks variable existence", {
  test_data <- create_minimal_test_data(
    n_subjects = 3, n_times = 3, with_covariates = FALSE
  )
  longitudinal_data <- test_data$longitudinal
  survival_data <- test_data$survival

  longitudinal_formula <- v ~ 1
  survival_formula <- Surv(time, status) ~ 1

  # Test missing ID variable
  expect_error(
    .validate(
      longitudinal_formula, longitudinal_data, survival_formula, survival_data,
      "subject_id", "time"
    ),
    "ID variable 'subject_id' not found in longitudinal data"
  )

  # Test missing time variable
  expect_error(
    .validate(
      longitudinal_formula, longitudinal_data, survival_formula, survival_data,
      "id", "timevar"
    ),
    "Time variable 'timevar' not found in longitudinal data"
  )

  # Test missing ID in survival data
  survival_data2 <- survival_data
  names(survival_data2)[1] <- "subject"
  expect_error(
    .validate(
      longitudinal_formula, longitudinal_data, survival_formula, survival_data2,
      "id", "time"
    ),
    "ID variable 'id' not found in survival data"
  )
})

test_that(".validate() checks formula variables", {
  test_data <- create_minimal_test_data(
    n_subjects = 3, n_times = 3, with_covariates = FALSE
  )
  longitudinal_data <- test_data$longitudinal
  survival_data <- test_data$survival

  # Test missing longitudinal formula variable
  expect_error(
    .validate(
      v ~ x1, longitudinal_data, Surv(time, status) ~ 1, survival_data,
      "id", "time"
    ),
    "Variables in longitudinal formula not found in data: x1"
  )

  # Test missing survival formula variable
  expect_error(
    .validate(
      v ~ 1, longitudinal_data, Surv(time, status) ~ w1, survival_data,
      "id", "time"
    ),
    "Variables in survival formula not found in data: w1"
  )
})

test_that(".validate() checks survival formula structure", {
  test_data <- create_minimal_test_data(
    n_subjects = 3, n_times = 3, with_covariates = FALSE
  )
  longitudinal_data <- test_data$longitudinal
  survival_data <- test_data$survival

  # Test invalid survival formula (no Surv())
  expect_error(
    .validate(v ~ 1, longitudinal_data, time ~ 1, survival_data, "id", "time"),
    "Survival formula must have Surv\\(\\) on the left-hand side"
  )

  # Test Surv() with missing variables
  expect_error(
    .validate(
      v ~ 1, longitudinal_data, Surv(time2, event) ~ 1, survival_data,
      "id", "time"
    ),
    "Variables in survival formula not found in data: time2, event"
  )
})

test_that(".validate() checks ID consistency", {
  # IDs 1-4 in longitudinal, only 1-3 in survival
  test_longitudinal_data <- create_minimal_test_data(
    n_subjects = 4, n_times = 3
  )
  test_survival_data <- create_minimal_test_data(n_subjects = 3, n_times = 3)
  longitudinal_data <- test_longitudinal_data$longitudinal
  survival_data <- test_survival_data$survival

  # Test orphaned longitudinal data
  expect_error(
    .validate(
      v ~ 1, longitudinal_data, Surv(time, status) ~ 1, survival_data,
      "id", "time"
    ),
    "Subjects in longitudinal data not found in survival data: 4"
  )

  # Test missing longitudinal data (should warn, not error)
  longitudinal_data2 <- data.frame(
    id = rep(1:2, each = 3), # Only IDs 1-2
    time = rep(0:2, 2),
    v = rnorm(6)
  )

  expect_warning(
    .validate(
      v ~ 1, longitudinal_data2, Surv(time, status) ~ 1, survival_data,
      "id", "time"
    ),
    "Subjects in survival data without longitudinal data: 3"
  )
})

test_that(".validate() checks time values", {
  # Test negative longitudinal times
  test_data <- create_minimal_test_data(n_subjects = 3, n_times = 3)
  longitudinal_data <- test_data$longitudinal
  longitudinal_data$time <- c(-1, 0, 1, 0, 1, 2, -0.5, 0.5, 1.5)
  survival_data <- test_data$survival

  expect_error(
    .validate(
      v ~ 1, longitudinal_data, Surv(time, status) ~ 1, survival_data,
      "id", "time"
    ),
    "Negative time values found in longitudinal data"
  )

  # Test non-positive survival times
  longitudinal_data2 <- data.frame(
    id = rep(1:3, each = 3),
    time = rep(0:2, 3),
    v = rnorm(9)
  )

  survival_data2 <- data.frame(
    id = 1:3,
    time = c(0, 3.0, 2.8), # 0 is not valid
    status = c(1, 0, 1)
  )

  expect_error(
    .validate(
      v ~ 1, longitudinal_data2, Surv(time, status) ~ 1, survival_data2,
      "id", "time"
    ),
    "Invalid observation times in survival data"
  )
})

test_that(".validate() checks status values", {
  test_data <- create_minimal_test_data(n_subjects = 3, n_times = 3)
  longitudinal_data <- test_data$longitudinal

  # Test invalid status values
  survival_data <- test_data$survival
  survival_data$status <- c(1, 2, 0) # 2 is invalid

  expect_error(
    .validate(
      v ~ 1, longitudinal_data, Surv(time, status) ~ 1, survival_data,
      "id", "time"
    ),
    "Invalid status values found: 2"
  )

  # Test missing status values
  survival_data2 <- data.frame(
    id = 1:3,
    time = c(2.5, 3.0, 2.8),
    status = c(1, NA, 0)
  )

  expect_error(
    .validate(
      v ~ 1, longitudinal_data, Surv(time, status) ~ 1, survival_data2,
      "id", "time"
    ),
    "Missing values in status variable"
  )
})

test_that(".validate() checks for NA values in critical columns", {
  # Test NA in longitudinal ID
  longitudinal_data <- data.frame(
    id = c(1, NA, 1, 2, 2, 2),
    time = c(0, 1, 2, 0, 1, 2),
    v = rnorm(6)
  )

  survival_data <- data.frame(
    id = 1:2,
    time = c(2.5, 3.0),
    status = c(1, 0)
  )

  expect_error(
    .validate(
      v ~ 1, longitudinal_data, Surv(time, status) ~ 1, survival_data,
      "id", "time"
    ),
    "Missing values found in ID in longitudinal data"
  )

  # Test NA in time variable
  longitudinal_data2 <- data.frame(
    id = rep(1:2, each = 3),
    time = c(0, NA, 2, 0, 1, 2),
    v = rnorm(6)
  )

  expect_error(
    .validate(
      v ~ 1, longitudinal_data2, Surv(time, status) ~ 1, survival_data,
      "id", "time"
    ),
    "Missing values found in Time in longitudinal data"
  )
})

test_that(".validate() checks for duplicate survival records", {
  longitudinal_data <- data.frame(
    id = rep(1:3, each = 3),
    time = rep(0:2, 3),
    v = rnorm(9)
  )

  survival_data <- data.frame(
    id = c(1, 2, 2, 3), # ID 2 is duplicated
    time = c(2.5, 3.0, 2.8, 3.2),
    status = c(1, 0, 1, 1)
  )

  expect_error(
    .validate(
      v ~ 1, longitudinal_data, Surv(time, status) ~ 1, survival_data,
      "id", "time"
    ),
    "Duplicate IDs found in survival data"
  )
})

test_that(".validate() checks data dimensions", {
  # Test empty longitudinal data
  longitudinal_data_empty <- data.frame(
    id = numeric(0),
    time = numeric(0),
    v = numeric(0)
  )

  survival_data <- data.frame(
    id = 1:3,
    time = c(2.5, 3.0, 2.8),
    status = c(1, 0, 1)
  )

  expect_error(
    .validate(
      v ~ 1, longitudinal_data_empty, Surv(time, status) ~ 1, survival_data,
      "id", "time"
    ),
    "Longitudinal data has no rows"
  )

  # Test empty survival data
  longitudinal_data <- data.frame(
    id = rep(1:3, each = 3),
    time = rep(0:2, 3),
    v = rnorm(9)
  )

  survival_data_empty <- data.frame(
    id = numeric(0),
    time = numeric(0),
    status = numeric(0)
  )

  expect_error(
    .validate(
      v ~ 1, longitudinal_data, Surv(time, status) ~ 1, survival_data_empty,
      "id", "time"
    ),
    "Survival data has no rows"
  )
})

test_that(".validate() warns about single observations per subject", {
  # Each subject has only one observation
  test_data <- create_minimal_test_data(n_subjects = 3, n_times = 1)
  longitudinal_data <- test_data$longitudinal
  survival_data <- test_data$survival

  expect_warning(
    .validate(
      v ~ 1, longitudinal_data, Surv(time, status) ~ 1, survival_data,
      "id", "time"
    ),
    "Each subject has only one longitudinal observation"
  )
})

test_that(".validate() works with valid data", {
  # Create completely valid data
  longitudinal_data <- create_test_longitudinal_data(
    n_subjects = 5, n_times = 4, seed = 123
  )[, c("id", "time", "v", "x1")]
  survival_data <- create_test_survival_data(n_subjects = 5, seed = 123)
  # Ensure times > max(long times)
  survival_data$time <- c(3.5, 4.0, 3.8, 4.2, 3.9)
  survival_data$w1 <- rnorm(5)

  # Should run without errors or warnings
  expect_silent(
    .validate(
      v ~ x1, longitudinal_data,  Surv(time, status) ~ w1, survival_data,
      "id", "time"
    )
  )
})

# ===========================================================================
# Additional tests for increased coverage
# ===========================================================================

test_that(".validate() handles more than 5 orphaned IDs in error message", {
  longitudinal_data_many <- data.frame(
    id = rep(1:10, each = 2),
    time = rep(0:1, 10),
    v = rnorm(20)
  )

  survival_data_few <- data.frame(
    id = 1:3, # Only 3 subjects in survival
    time = c(2, 2.5, 3),
    status = c(1, 0, 1)
  )

  expect_error(
    .validate(
      v ~ 1, longitudinal_data_many, Surv(time, status) ~ 1, survival_data_few,
      "id", "time"
    ),
    "Subjects in longitudinal data not found in survival data: 4, 5, 6, 7, 8"
  )
})

test_that(".validate() warns about more than 5 missing longitudinal subjects", {
  longitudinal_data_few <- data.frame(
    id = rep(1:2, each = 2),
    time = rep(0:1, 2),
    v = rnorm(4)
  )

  survival_data_many <- data.frame(
    id = 1:10, # 10 subjects in survival
    time = rexp(10, 0.1) + 0.1,
    status = rbinom(10, 1, 0.7)
  )

  expect_warning(
    .validate(
      v ~ 1, longitudinal_data_few, Surv(time, status) ~ 1, survival_data_many,
      "id", "time"
    ),
    "Subjects in survival data without longitudinal data: 3, 4, 5, 6, 7"
  )
})

test_that(".validate() handles boundary conditions", {
  # Minimum viable data (2 subjects, 2 time points each)
  longitudinal_data_min <- data.frame(
    id = c(1, 1, 2, 2),
    time = c(0, 1, 0, 1),
    v = c(1, 2, 3, 4)
  )

  survival_data_min <- data.frame(
    id = c(1, 2),
    time = c(2, 2.5),
    status = c(1, 0)
  )

  expect_silent(
    .validate(
      v ~ 1, longitudinal_data_min, Surv(time, status) ~ 1, survival_data_min,
      "id", "time"
    )
  )
})

test_that(".validate() handles special characters in variable names", {
  # Data with special characters in column names
  longitudinal_data_special <- data.frame(
    `subject.id` = rep(1:3, each = 3),
    `time_point` = rep(0:2, 3),
    `outcome-1` = rnorm(9),
    check.names = FALSE
  )

  survival_data_special <- data.frame(
    `subject.id` = 1:3,
    `event.time` = c(2.5, 3.0, 2.8),
    status = c(1, 0, 1),
    check.names = FALSE
  )

  expect_silent(
    .validate(
      `outcome-1` ~ 1, longitudinal_data_special,
      Surv(`event.time`, status) ~ 1, survival_data_special,
      "subject.id", "time_point"
    )
  )
})

test_that(".validate() handles large datasets efficiently", {
  set.seed(42)
  n_subjects <- 1000
  n_times <- 10

  # Large longitudinal dataset
  longitudinal_data_large <- data.frame(
    id = rep(1:n_subjects, each = n_times),
    time = rep(0:(n_times - 1), n_subjects),
    v = rnorm(n_subjects * n_times),
    x1 = rnorm(n_subjects * n_times),
    x2 = rnorm(n_subjects * n_times)
  )

  # Large survival dataset
  survival_data_large <- data.frame(
    id = 1:n_subjects,
    time = rexp(n_subjects, rate = 0.1) + 0.1,
    status = rbinom(n_subjects, 1, 0.7),
    w1 = rnorm(n_subjects),
    w2 = rbinom(n_subjects, 1, 0.5)
  )

  # Should complete within reasonable time
  time_taken <- system.time({
    suppressWarnings(
      .validate(
        v ~ x1 + x2,  longitudinal_data_large,
        Surv(time, status) ~ w1 + w2, survival_data_large,
        "id", "time"
      )
    )
  })

  expect_lt(time_taken["elapsed"], 5)
})

test_that(".validate() checks spline_baseline parameters", {
  test_data <- create_minimal_test_data(n_subjects = 3, n_times = 3)
  longitudinal_data <- test_data$longitudinal
  survival_data <- test_data$survival

  # Test invalid parameter names
  expect_error(
    .validate(
      v ~ 1, longitudinal_data, Surv(time, status) ~ 1,  survival_data,
      "id", "time", spline_baseline = list(invalid_param = 3),
      spline_index = list()
    ),
    "Invalid parameters in spline_baseline: invalid_param"
  )

  # Test invalid degree values
  expect_error(
    .validate(
      v ~ 1, longitudinal_data, Surv(time, status) ~ 1, survival_data,
      "id", "time", spline_baseline = list(degree = 0), spline_index = list()
    ),
    "spline_baseline\\$degree must be a single integer between 1 and 5"
  )

  expect_error(
    .validate(
      v ~ 1, longitudinal_data, Surv(time, status) ~ 1, survival_data,
      "id", "time", spline_baseline = list(degree = 6), spline_index = list()
    ),
    "spline_baseline\\$degree must be a single integer between 1 and 5"
  )

  # Test invalid n_knots values
  expect_error(
    .validate(
      v ~ 1, longitudinal_data, Surv(time, status) ~ 1, survival_data,
      "id", "time", spline_baseline = list(n_knots = 0), spline_index = list()
    ),
    "spline_baseline\\$n_knots must be a single integer between 1 and 20"
  )

  expect_error(
    .validate(
      v ~ 1, longitudinal_data, Surv(time, status) ~ 1, survival_data,
      "id", "time", spline_baseline = list(n_knots = 21), spline_index = list()
    ),
    "spline_baseline\\$n_knots must be a single integer between 1 and 20"
  )

  # Test invalid knot_placement values
  expect_error(
    .validate(
      v ~ 1, longitudinal_data, Surv(time, status) ~ 1, survival_data,
      "id", "time", spline_baseline = list(knot_placement = "invalid"),
      spline_index = list()
    ),
    "spline_baseline\\$knot_placement must be one of: quantile, equal"
  )

  # Test invalid boundary_knots
  expect_error(
    .validate(
      v ~ 1, longitudinal_data, Surv(time, status) ~ 1, survival_data,
      "id", "time", spline_baseline = list(boundary_knots = c(1)),
      spline_index = list()
    ),
    paste(
      "spline_baseline\\$boundary_knots must be NULL or",
      "a numeric vector of length 2"
    )
  )

  expect_error(
    .validate(
      v ~ 1, longitudinal_data, Surv(time, status) ~ 1, survival_data,
      "id", "time", spline_baseline = list(boundary_knots = c(2, 1)),
      spline_index = list()
    ),
    paste(
      "spline_baseline\\$boundary_knots\\[1\\] must be",
      "less than boundary_knots\\[2\\]"
    )
  )

  # Test valid spline_baseline parameters
  expect_silent(
    .validate(
      v ~ 1, longitudinal_data, Surv(time, status) ~ 1, survival_data,
      "id", "time",
      spline_baseline = list(
        degree = 3,
        n_knots = 5,
        knot_placement = "quantile",
        boundary_knots = c(0, 10)
      ),
      spline_index = list()
    )
  )
})

test_that(".validate() checks spline_index parameters", {
  test_data <- create_minimal_test_data(n_subjects = 3, n_times = 3)
  longitudinal_data <- test_data$longitudinal
  survival_data <- test_data$survival

  # Test invalid parameter names
  expect_error(
    .validate(
      v ~ 1, longitudinal_data, Surv(time, status) ~ 1, survival_data,
      "id", "time", spline_baseline = list(),
      spline_index = list(invalid_param = 3)
    ),
    "Invalid parameters in spline_index: invalid_param"
  )

  # Test invalid degree values
  expect_error(
    .validate(
      v ~ 1, longitudinal_data, Surv(time, status) ~ 1, survival_data,
      "id", "time", spline_baseline = list(), spline_index = list(degree = 0)
    ),
    "spline_index\\$degree must be a single integer between 1 and 5"
  )

  expect_error(
    .validate(
      v ~ 1, longitudinal_data, Surv(time, status) ~ 1, survival_data,
      "id", "time", spline_baseline = list(), spline_index = list(degree = 6)
    ),
    "spline_index\\$degree must be a single integer between 1 and 5"
  )

  # Test invalid n_knots values
  expect_error(
    .validate(
      v ~ 1, longitudinal_data, Surv(time, status) ~ 1, survival_data,
      "id", "time", spline_baseline = list(), spline_index = list(n_knots = 0)
    ),
    "spline_index\\$n_knots must be a single integer between 1 and 20"
  )

  expect_error(
    .validate(
      v ~ 1, longitudinal_data, Surv(time, status) ~ 1, survival_data,
      "id", "time", spline_baseline = list(), spline_index = list(n_knots = 21)
    ),
    "spline_index\\$n_knots must be a single integer between 1 and 20"
  )

  # Test invalid knot_placement values
  expect_error(
    .validate(
      v ~ 1, longitudinal_data, Surv(time, status) ~ 1, survival_data,
      "id", "time", spline_baseline = list(),
      spline_index = list(knot_placement = "invalid")
    ),
    "spline_index\\$knot_placement must be one of: quantile, equal"
  )

  # Test invalid boundary_knots
  expect_error(
    .validate(
      v ~ 1, longitudinal_data, Surv(time, status) ~ 1, survival_data,
      "id", "time", spline_baseline = list(),
      spline_index = list(boundary_knots = c(1, 2, 3))
    ),
    "spline_index\\$boundary_knots must be NULL or a numeric vector of length 2"
  )

  expect_error(
    .validate(
      v ~ 1, longitudinal_data, Surv(time, status) ~ 1,  survival_data,
      "id", "time", spline_baseline = list(),
      spline_index = list(boundary_knots = c(5, 5))
    ),
    paste(
      "spline_index\\$boundary_knots\\[1\\] must be",
      "less than boundary_knots\\[2\\]"
    )
  )

  # Test valid spline_index parameters
  expect_silent(
    .validate(
      v ~ 1, longitudinal_data, Surv(time, status) ~ 1, survival_data,
      "id", "time", spline_baseline = list(),
      spline_index = list(
        degree = 2,
        n_knots = 4,
        knot_placement = "equal",
        boundary_knots = c(-1, 5)
      )
    )
  )
})

test_that(".validate() checks that spline parameters are lists", {
  test_data <- create_minimal_test_data(n_subjects = 3, n_times = 3)
  longitudinal_data <- test_data$longitudinal
  survival_data <- test_data$survival

  # Test non-list spline_baseline
  expect_error(
    .validate(
      v ~ 1, longitudinal_data, Surv(time, status) ~ 1, survival_data,
      "id", "time", spline_baseline = "not a list", spline_index = list()
    ),
    "spline_baseline must be a list"
  )

  # Test non-list spline_index
  expect_error(
    .validate(
      v ~ 1, longitudinal_data, Surv(time, status) ~ 1, survival_data,
      "id", "time", spline_baseline = list(), spline_index = 5
    ),
    "spline_index must be a list"
  )
})

test_that(".validate() handles different data types", {
  # Integer IDs
  longitudinal_data_int <- data.frame(
    id = as.integer(rep(1:3, each = 3)),
    time = rep(0:2, 3),
    v = rnorm(9)
  )

  survival_data_int <- data.frame(
    id = as.integer(1:3),
    time = c(2.5, 3.0, 2.8),
    status = c(1, 0, 1)
  )

  expect_silent(
    .validate(
      v ~ 1,
      longitudinal_data_int,
      Surv(time, status) ~ 1,
      survival_data_int,
      "id",
      "time"
    )
  )

  # Character IDs
  longitudinal_data_char <- data.frame(
    id = rep(c("A", "B", "C"), each = 3),
    time = rep(0:2, 3),
    v = rnorm(9)
  )

  survival_data_char <- data.frame(
    id = c("A", "B", "C"),
    time = c(2.5, 3.0, 2.8),
    status = c(1, 0, 1)
  )

  expect_silent(
    .validate(
      v ~ 1,
      longitudinal_data_char,
      Surv(time, status) ~ 1,
      survival_data_char,
      "id",
      "time"
    )
  )
})
