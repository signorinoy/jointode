test_that("validate() checks input types correctly", {
  # Create valid test data
  data.long <- create_test_longitudinal_data(
    n_subjects = 3, n_times = 3, seed = 123
  )[, c("id", "time", "v", "x1")]
  data.surv <- create_test_survival_data(
    n_subjects = 3, seed = 123
  )[, c("id", "time", "status")]
  data.surv$w1 <- rnorm(3)

  formula.long <- v ~ x1
  formula.surv <- Surv(time, status) ~ w1

  # Test invalid formula types
  expect_error(
    validate("not a formula", formula.surv, data.long, data.surv, "id", "time"),
    "formula.long must be a formula"
  )

  expect_error(
    validate(formula.long, "not a formula", data.long, data.surv, "id", "time"),
    "formula.surv must be a formula"
  )

  # Test invalid data frame types
  expect_error(
    validate(formula.long, formula.surv, "not a df", data.surv, "id", "time"),
    "data.long must be a data frame"
  )

  expect_error(
    validate(formula.long, formula.surv, data.long, list(), "id", "time"),
    "data.surv must be a data frame"
  )

  # Test invalid ID and time types
  expect_error(
    validate(formula.long, formula.surv, data.long, data.surv, 123, "time"),
    "id must be a single character string"
  )

  expect_error(
    validate(
      formula.long, formula.surv, data.long, data.surv,
      "id", c("t1", "t2")
    ),
    "time must be a single character string"
  )
})

test_that("validate() checks variable existence", {
  test_data <- create_minimal_test_data(
    n_subjects = 3, n_times = 3, with_covariates = FALSE
  )
  data.long <- test_data$longitudinal
  data.surv <- test_data$survival

  formula.long <- v ~ 1
  formula.surv <- Surv(time, status) ~ 1

  # Test missing ID variable
  expect_error(
    validate(
      formula.long, formula.surv, data.long, data.surv,
      "subject_id", "time"
    ),
    "ID variable 'subject_id' not found in longitudinal data"
  )

  # Test missing time variable
  expect_error(
    validate(formula.long, formula.surv, data.long, data.surv, "id", "timevar"),
    "Time variable 'timevar' not found in longitudinal data"
  )

  # Test missing ID in survival data
  data.surv2 <- data.surv
  names(data.surv2)[1] <- "subject"
  expect_error(
    validate(formula.long, formula.surv, data.long, data.surv2, "id", "time"),
    "ID variable 'id' not found in survival data"
  )
})

test_that("validate() checks formula variables", {
  test_data <- create_minimal_test_data(
    n_subjects = 3, n_times = 3, with_covariates = FALSE
  )
  data.long <- test_data$longitudinal
  data.surv <- test_data$survival

  # Test missing longitudinal formula variable
  expect_error(
    validate(
      v ~ x1, Surv(time, status) ~ 1,
      data.long, data.surv, "id", "time"
    ),
    "Variables in longitudinal formula not found in data: x1"
  )

  # Test missing survival formula variable
  expect_error(
    validate(
      v ~ 1, Surv(time, status) ~ w1,
      data.long, data.surv, "id", "time"
    ),
    "Variables in survival formula not found in data: w1"
  )
})

test_that("validate() checks survival formula structure", {
  test_data <- create_minimal_test_data(
    n_subjects = 3, n_times = 3, with_covariates = FALSE
  )
  data.long <- test_data$longitudinal
  data.surv <- test_data$survival

  # Test invalid survival formula (no Surv())
  expect_error(
    validate(v ~ 1, time ~ 1, data.long, data.surv, "id", "time"),
    "Survival formula must have Surv\\(\\) on the left-hand side"
  )

  # Test Surv() with missing variables
  expect_error(
    validate(
      v ~ 1, Surv(time2, event) ~ 1, data.long, data.surv,
      "id", "time"
    ),
    "Variables in survival formula not found in data: time2, event"
  )
})

test_that("validate() checks ID consistency", {
  # IDs 1-4 in longitudinal, only 1-3 in survival
  test_data_long <- create_minimal_test_data(n_subjects = 4, n_times = 3)
  test_data_surv <- create_minimal_test_data(n_subjects = 3, n_times = 3)
  data.long <- test_data_long$longitudinal
  data.surv <- test_data_surv$survival

  # Test orphaned longitudinal data
  expect_error(
    validate(
      v ~ 1, Surv(time, status) ~ 1,
      data.long, data.surv, "id", "time"
    ),
    "Subjects in longitudinal data not found in survival data: 4"
  )

  # Test missing longitudinal data (should warn, not error)
  data.long2 <- data.frame(
    id = rep(1:2, each = 3), # Only IDs 1-2
    time = rep(0:2, 2),
    v = rnorm(6)
  )

  expect_warning(
    validate(
      v ~ 1, Surv(time, status) ~ 1,
      data.long2, data.surv, "id", "time"
    ),
    "Subjects in survival data without longitudinal data: 3"
  )
})

test_that("validate() checks time values", {
  # Test negative longitudinal times
  test_data <- create_minimal_test_data(n_subjects = 3, n_times = 3)
  data.long <- test_data$longitudinal
  data.long$time <- c(-1, 0, 1, 0, 1, 2, -0.5, 0.5, 1.5)
  data.surv <- test_data$survival

  expect_error(
    validate(
      v ~ 1, Surv(time, status) ~ 1,
      data.long, data.surv, "id", "time"
    ),
    "Negative time values found in longitudinal data"
  )

  # Test non-positive survival times
  data.long2 <- data.frame(
    id = rep(1:3, each = 3),
    time = rep(0:2, 3),
    v = rnorm(9)
  )

  data.surv2 <- data.frame(
    id = 1:3,
    time = c(0, 3.0, 2.8), # 0 is not valid
    status = c(1, 0, 1)
  )

  expect_error(
    validate(
      v ~ 1, Surv(time, status) ~ 1,
      data.long2, data.surv2, "id", "time"
    ),
    "Invalid observation times in survival data"
  )
})

test_that("validate() checks status values", {
  test_data <- create_minimal_test_data(n_subjects = 3, n_times = 3)
  data.long <- test_data$longitudinal

  # Test invalid status values
  data.surv <- test_data$survival
  data.surv$status <- c(1, 2, 0) # 2 is invalid

  expect_error(
    validate(
      v ~ 1, Surv(time, status) ~ 1,
      data.long, data.surv, "id", "time"
    ),
    "Invalid status values found: 2"
  )

  # Test missing status values
  data.surv2 <- data.frame(
    id = 1:3,
    time = c(2.5, 3.0, 2.8),
    status = c(1, NA, 0)
  )

  expect_error(
    validate(
      v ~ 1, Surv(time, status) ~ 1,
      data.long, data.surv2, "id", "time"
    ),
    "Missing values in status variable"
  )
})

test_that("validate() checks for NA values in critical columns", {
  # Test NA in longitudinal ID
  data.long <- data.frame(
    id = c(1, NA, 1, 2, 2, 2),
    time = c(0, 1, 2, 0, 1, 2),
    v = rnorm(6)
  )

  data.surv <- data.frame(
    id = 1:2,
    time = c(2.5, 3.0),
    status = c(1, 0)
  )

  expect_error(
    validate(
      v ~ 1, Surv(time, status) ~ 1,
      data.long, data.surv, "id", "time"
    ),
    "Missing values found in ID in longitudinal data"
  )

  # Test NA in time variable
  data.long2 <- data.frame(
    id = rep(1:2, each = 3),
    time = c(0, NA, 2, 0, 1, 2),
    v = rnorm(6)
  )

  expect_error(
    validate(
      v ~ 1, Surv(time, status) ~ 1,
      data.long2, data.surv, "id", "time"
    ),
    "Missing values found in Time in longitudinal data"
  )
})

test_that("validate() checks for duplicate survival records", {
  data.long <- data.frame(
    id = rep(1:3, each = 3),
    time = rep(0:2, 3),
    v = rnorm(9)
  )

  data.surv <- data.frame(
    id = c(1, 2, 2, 3), # ID 2 is duplicated
    time = c(2.5, 3.0, 2.8, 3.2),
    status = c(1, 0, 1, 1)
  )

  expect_error(
    validate(
      v ~ 1, Surv(time, status) ~ 1,
      data.long, data.surv, "id", "time"
    ),
    "Duplicate IDs found in survival data"
  )
})

test_that("validate() checks data dimensions", {
  # Test empty longitudinal data
  data_long_empty <- data.frame(
    id = numeric(0),
    time = numeric(0),
    v = numeric(0)
  )

  data.surv <- data.frame(
    id = 1:3,
    time = c(2.5, 3.0, 2.8),
    status = c(1, 0, 1)
  )

  expect_error(
    validate(
      v ~ 1, Surv(time, status) ~ 1,
      data_long_empty, data.surv, "id", "time"
    ),
    "Longitudinal data has no rows"
  )

  # Test empty survival data
  data.long <- data.frame(
    id = rep(1:3, each = 3),
    time = rep(0:2, 3),
    v = rnorm(9)
  )

  data_surv_empty <- data.frame(
    id = numeric(0),
    time = numeric(0),
    status = numeric(0)
  )

  expect_error(
    validate(
      v ~ 1, Surv(time, status) ~ 1,
      data.long, data_surv_empty, "id", "time"
    ),
    "Survival data has no rows"
  )
})

test_that("validate() warns about single observations per subject", {
  # Each subject has only one observation
  test_data <- create_minimal_test_data(n_subjects = 3, n_times = 1)
  data.long <- test_data$longitudinal
  data.surv <- test_data$survival

  expect_warning(
    validate(
      v ~ 1, Surv(time, status) ~ 1,
      data.long, data.surv, "id", "time"
    ),
    "Each subject has only one longitudinal observation"
  )
})

test_that("validate() works with valid data", {
  # Create completely valid data
  data.long <- create_test_longitudinal_data(
    n_subjects = 5, n_times = 4, seed = 123
  )[, c("id", "time", "v", "x1")]
  data.surv <- create_test_survival_data(n_subjects = 5, seed = 123)
  data.surv$time <- c(3.5, 4.0, 3.8, 4.2, 3.9) # Ensure times > max(long times)
  data.surv$w1 <- rnorm(5)

  # Should run without errors or warnings
  expect_silent(
    validate(
      v ~ x1, Surv(time, status) ~ w1,
      data.long, data.surv, "id", "time"
    )
  )
})

# ===========================================================================
# Additional tests for increased coverage
# ===========================================================================

test_that("validate() handles more than 5 orphaned IDs in error message", {
  data.long.many <- data.frame(
    id = rep(1:10, each = 2),
    time = rep(0:1, 10),
    v = rnorm(20)
  )

  data.surv.few <- data.frame(
    id = 1:3, # Only 3 subjects in survival
    time = c(2, 2.5, 3),
    status = c(1, 0, 1)
  )

  expect_error(
    validate(
      v ~ 1,
      Surv(time, status) ~ 1,
      data.long.many,
      data.surv.few,
      "id",
      "time"
    ),
    "Subjects in longitudinal data not found in survival data: 4, 5, 6, 7, 8"
  )
})

test_that("validate() warns about more than 5 missing longitudinal subjects", {
  data.long.few <- data.frame(
    id = rep(1:2, each = 2),
    time = rep(0:1, 2),
    v = rnorm(4)
  )

  data.surv.many <- data.frame(
    id = 1:10, # 10 subjects in survival
    time = rexp(10, 0.1) + 0.1,
    status = rbinom(10, 1, 0.7)
  )

  expect_warning(
    validate(
      v ~ 1,
      Surv(time, status) ~ 1,
      data.long.few,
      data.surv.many,
      "id",
      "time"
    ),
    "Subjects in survival data without longitudinal data: 3, 4, 5, 6, 7"
  )
})

test_that("validate() handles boundary conditions", {
  # Minimum viable data (2 subjects, 2 time points each)
  data.long.min <- data.frame(
    id = c(1, 1, 2, 2),
    time = c(0, 1, 0, 1),
    v = c(1, 2, 3, 4)
  )

  data.surv.min <- data.frame(
    id = c(1, 2),
    time = c(2, 2.5),
    status = c(1, 0)
  )

  expect_silent(
    validate(
      v ~ 1,
      Surv(time, status) ~ 1,
      data.long.min,
      data.surv.min,
      "id",
      "time"
    )
  )
})

test_that("validate() handles special characters in variable names", {
  # Data with special characters in column names
  data.long.special <- data.frame(
    `subject.id` = rep(1:3, each = 3),
    `time_point` = rep(0:2, 3),
    `outcome-1` = rnorm(9),
    check.names = FALSE
  )

  data.surv.special <- data.frame(
    `subject.id` = 1:3,
    `event.time` = c(2.5, 3.0, 2.8),
    status = c(1, 0, 1),
    check.names = FALSE
  )

  expect_silent(
    validate(
      `outcome-1` ~ 1,
      Surv(`event.time`, status) ~ 1,
      data.long.special,
      data.surv.special,
      "subject.id",
      "time_point"
    )
  )
})

test_that("validate() handles large datasets efficiently", {
  set.seed(42)
  n_subjects <- 1000
  n_times <- 10

  # Large longitudinal dataset
  data.long.large <- data.frame(
    id = rep(1:n_subjects, each = n_times),
    time = rep(0:(n_times - 1), n_subjects),
    v = rnorm(n_subjects * n_times),
    x1 = rnorm(n_subjects * n_times),
    x2 = rnorm(n_subjects * n_times)
  )

  # Large survival dataset
  data.surv.large <- data.frame(
    id = 1:n_subjects,
    time = rexp(n_subjects, rate = 0.1) + 0.1,
    status = rbinom(n_subjects, 1, 0.7),
    w1 = rnorm(n_subjects),
    w2 = rbinom(n_subjects, 1, 0.5)
  )

  # Should complete within reasonable time
  time_taken <- system.time({
    suppressWarnings(
      validate(
        v ~ x1 + x2,
        Surv(time, status) ~ w1 + w2,
        data.long.large,
        data.surv.large,
        "id",
        "time"
      )
    )
  })

  expect_lt(time_taken["elapsed"], 5)
})

test_that("validate() handles different data types", {
  # Integer IDs
  data.long.int <- data.frame(
    id = as.integer(rep(1:3, each = 3)),
    time = rep(0:2, 3),
    v = rnorm(9)
  )

  data.surv.int <- data.frame(
    id = as.integer(1:3),
    time = c(2.5, 3.0, 2.8),
    status = c(1, 0, 1)
  )

  expect_silent(
    validate(
      v ~ 1,
      Surv(time, status) ~ 1,
      data.long.int,
      data.surv.int,
      "id",
      "time"
    )
  )

  # Character IDs
  data.long.char <- data.frame(
    id = rep(c("A", "B", "C"), each = 3),
    time = rep(0:2, 3),
    v = rnorm(9)
  )

  data.surv.char <- data.frame(
    id = c("A", "B", "C"),
    time = c(2.5, 3.0, 2.8),
    status = c(1, 0, 1)
  )

  expect_silent(
    validate(
      v ~ 1,
      Surv(time, status) ~ 1,
      data.long.char,
      data.surv.char,
      "id",
      "time"
    )
  )
})
