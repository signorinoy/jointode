# Load required packages
library(survival)

test_that("validate() checks input types correctly", {
  # Create valid test data
  set.seed(123)
  data.long <- data.frame(
    id = rep(1:3, each = 3),
    time = rep(0:2, 3),
    y = rnorm(9),
    x1 = rnorm(9)
  )

  data.surv <- data.frame(
    id = 1:3,
    time = c(2.5, 3.0, 2.8),
    status = c(1, 0, 1),
    x2 = rnorm(3)
  )

  formula.long <- y ~ time + x1
  formula.surv <- Surv(time, status) ~ x2

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
  data.long <- data.frame(
    id = rep(1:3, each = 3),
    time = rep(0:2, 3),
    y = rnorm(9)
  )

  data.surv <- data.frame(
    id = 1:3,
    time = c(2.5, 3.0, 2.8),
    status = c(1, 0, 1)
  )

  formula.long <- y ~ time
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
  data.long <- data.frame(
    id = rep(1:3, each = 3),
    time = rep(0:2, 3),
    y = rnorm(9)
  )

  data.surv <- data.frame(
    id = 1:3,
    time = c(2.5, 3.0, 2.8),
    status = c(1, 0, 1)
  )

  # Test missing longitudinal formula variable
  expect_error(
    validate(
      y ~ time + x1, Surv(time, status) ~ 1,
      data.long, data.surv, "id", "time"
    ),
    "Variables in longitudinal formula not found in data: x1"
  )

  # Test missing survival formula variable
  expect_error(
    validate(
      y ~ time, Surv(time, status) ~ x2,
      data.long, data.surv, "id", "time"
    ),
    "Variables in survival formula not found in data: x2"
  )
})

test_that("validate() checks survival formula structure", {
  data.long <- data.frame(
    id = rep(1:3, each = 3),
    time = rep(0:2, 3),
    y = rnorm(9)
  )

  data.surv <- data.frame(
    id = 1:3,
    time = c(2.5, 3.0, 2.8),
    status = c(1, 0, 1)
  )

  # Test invalid survival formula (no Surv())
  expect_error(
    validate(y ~ time, time ~ 1, data.long, data.surv, "id", "time"),
    "Survival formula must have Surv\\(\\) on the left-hand side"
  )

  # Test Surv() with missing variables
  expect_error(
    validate(
      y ~ time, Surv(time2, event) ~ 1, data.long, data.surv,
      "id", "time"
    ),
    "Variables in survival formula not found in data: time2, event"
  )
})

test_that("validate() checks ID consistency", {
  data.long <- data.frame(
    id = rep(1:4, each = 3), # IDs 1-4
    time = rep(0:2, 4),
    y = rnorm(12)
  )

  data.surv <- data.frame(
    id = 1:3, # Only IDs 1-3
    time = c(2.5, 3.0, 2.8),
    status = c(1, 0, 1)
  )

  # Test orphaned longitudinal data
  expect_error(
    validate(
      y ~ time, Surv(time, status) ~ 1,
      data.long, data.surv, "id", "time"
    ),
    "Subjects in longitudinal data not found in survival data: 4"
  )

  # Test missing longitudinal data (should warn, not error)
  data.long2 <- data.frame(
    id = rep(1:2, each = 3), # Only IDs 1-2
    time = rep(0:2, 2),
    y = rnorm(6)
  )

  expect_warning(
    validate(
      y ~ time, Surv(time, status) ~ 1,
      data.long2, data.surv, "id", "time"
    ),
    "Subjects in survival data without longitudinal data: 3"
  )
})

test_that("validate() checks time values", {
  # Test negative longitudinal times
  data.long <- data.frame(
    id = rep(1:3, each = 3),
    time = c(-1, 0, 1, 0, 1, 2, -0.5, 0.5, 1.5),
    y = rnorm(9)
  )

  data.surv <- data.frame(
    id = 1:3,
    time = c(2.5, 3.0, 2.8),
    status = c(1, 0, 1)
  )

  expect_error(
    validate(
      y ~ time, Surv(time, status) ~ 1,
      data.long, data.surv, "id", "time"
    ),
    "Negative time values found in longitudinal data"
  )

  # Test non-positive survival times
  data.long2 <- data.frame(
    id = rep(1:3, each = 3),
    time = rep(0:2, 3),
    y = rnorm(9)
  )

  data.surv2 <- data.frame(
    id = 1:3,
    time = c(0, 3.0, 2.8), # 0 is not valid
    status = c(1, 0, 1)
  )

  expect_error(
    validate(
      y ~ time, Surv(time, status) ~ 1,
      data.long2, data.surv2, "id", "time"
    ),
    "Invalid observation times in survival data"
  )
})

test_that("validate() checks status values", {
  data.long <- data.frame(
    id = rep(1:3, each = 3),
    time = rep(0:2, 3),
    y = rnorm(9)
  )

  # Test invalid status values
  data.surv <- data.frame(
    id = 1:3,
    time = c(2.5, 3.0, 2.8),
    status = c(1, 2, 0) # 2 is invalid
  )

  expect_error(
    validate(
      y ~ time, Surv(time, status) ~ 1,
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
      y ~ time, Surv(time, status) ~ 1,
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
    y = rnorm(6)
  )

  data.surv <- data.frame(
    id = 1:2,
    time = c(2.5, 3.0),
    status = c(1, 0)
  )

  expect_error(
    validate(
      y ~ time, Surv(time, status) ~ 1,
      data.long, data.surv, "id", "time"
    ),
    "Missing values found in ID in longitudinal data"
  )

  # Test NA in time variable
  data.long2 <- data.frame(
    id = rep(1:2, each = 3),
    time = c(0, NA, 2, 0, 1, 2),
    y = rnorm(6)
  )

  expect_error(
    validate(
      y ~ time, Surv(time, status) ~ 1,
      data.long2, data.surv, "id", "time"
    ),
    "Missing values found in Time in longitudinal data"
  )
})

test_that("validate() checks for duplicate survival records", {
  data.long <- data.frame(
    id = rep(1:3, each = 3),
    time = rep(0:2, 3),
    y = rnorm(9)
  )

  data.surv <- data.frame(
    id = c(1, 2, 2, 3), # ID 2 is duplicated
    time = c(2.5, 3.0, 2.8, 3.2),
    status = c(1, 0, 1, 1)
  )

  expect_error(
    validate(
      y ~ time, Surv(time, status) ~ 1,
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
    y = numeric(0)
  )

  data.surv <- data.frame(
    id = 1:3,
    time = c(2.5, 3.0, 2.8),
    status = c(1, 0, 1)
  )

  expect_error(
    validate(
      y ~ time, Surv(time, status) ~ 1,
      data_long_empty, data.surv, "id", "time"
    ),
    "Longitudinal data has no rows"
  )

  # Test empty survival data
  data.long <- data.frame(
    id = rep(1:3, each = 3),
    time = rep(0:2, 3),
    y = rnorm(9)
  )

  data_surv_empty <- data.frame(
    id = numeric(0),
    time = numeric(0),
    status = numeric(0)
  )

  expect_error(
    validate(
      y ~ time, Surv(time, status) ~ 1,
      data.long, data_surv_empty, "id", "time"
    ),
    "Survival data has no rows"
  )
})

test_that("validate() warns about single observations per subject", {
  data.long <- data.frame(
    id = 1:3, # Each subject has only one observation
    time = c(0, 0, 0),
    y = rnorm(3)
  )

  data.surv <- data.frame(
    id = 1:3,
    time = c(2.5, 3.0, 2.8),
    status = c(1, 0, 1)
  )

  expect_warning(
    validate(
      y ~ time, Surv(time, status) ~ 1,
      data.long, data.surv, "id", "time"
    ),
    "Each subject has only one longitudinal observation"
  )
})

test_that("validate() works with valid data", {
  # Create completely valid data
  set.seed(123)
  data.long <- data.frame(
    id = rep(1:5, each = 4),
    time = rep(0:3, 5),
    y = rnorm(20),
    x1 = rnorm(20)
  )

  data.surv <- data.frame(
    id = 1:5,
    time = c(3.5, 4.0, 3.8, 4.2, 3.9),
    status = c(1, 0, 1, 1, 0),
    x2 = rnorm(5)
  )

  # Should run without errors or warnings
  expect_silent(
    validate(
      y ~ time + x1, Surv(time, status) ~ x2,
      data.long, data.surv, "id", "time"
    )
  )
})
