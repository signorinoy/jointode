# Test process function
library(survival)

test_that("process() correctly structures basic data", {
  data.long <- create_test_longitudinal_data(
    n_subjects = 3, n_times = 4, seed = 123
  )
  data.long$v <- data.long$v + 10 # Adjust mean

  data.surv <- create_test_survival_data(n_subjects = 3, seed = 123)
  names(data.surv)[2] <- "obstime" # Rename time column

  result <- process(
    formula.long = v ~ x1 + x2,
    formula.surv = Surv(obstime, status) ~ w1 + w2,
    data.long = data.long,
    data.surv = data.surv,
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
  expect_true(attr(result, "event_rate") >= 0 &&
                attr(result, "event_rate") <= 1)
})

test_that("process() handles missing longitudinal data", {
  test_data <- create_mismatched_test_data(
    n_subjects_long = 2, n_subjects_surv = 3
  )
  data.long <- test_data$longitudinal
  data.surv <- test_data$survival
  names(data.surv)[2] <- "obstime"

  result <- process(
    formula.long = v ~ 1,
    formula.surv = Surv(obstime, status) ~ 1,
    data.long = data.long,
    data.surv = data.surv,
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

test_that("process() handles formula transformations", {
  set.seed(456)
  data.long <- data.frame(
    id = rep(1:10, each = 3),
    time = rep(1:3, 10),
    v = exp(rnorm(30)),
    x = runif(30, 1, 10)
  )

  data.surv <- data.frame(
    id = 1:10,
    surv_time = rexp(10),
    status = rbinom(10, 1, 0.5),
    z = runif(10, 1, 10)
  )

  result <- process(
    formula.long = log(v) ~ log(x) + I(x^2),
    formula.surv = Surv(surv_time, status) ~ log(z),
    data.long = data.long,
    data.surv = data.surv,
    id = "id",
    time = "time"
  )

  # Check transformations applied
  expect_equal(
    unname(result[[1]]$longitudinal$measurements[1]),
    log(data.long$v[data.long$id == 1][1])
  )
  expect_equal(ncol(result[[1]]$longitudinal$covariates), 2)
})

test_that("process() handles interaction terms", {
  data.long <- data.frame(
    id = rep(1:5, each = 3),
    time = rep(0:2, 5),
    v = rnorm(15),
    x1 = rnorm(15),
    x2 = rnorm(15)
  )

  data.surv <- data.frame(
    id = 1:5,
    surv_time = rexp(5),
    status = rbinom(5, 1, 0.5),
    w1 = rnorm(5),
    w2 = rnorm(5)
  )

  result <- process(
    formula.long = v ~ x1 * x2,
    formula.surv = Surv(surv_time, status) ~ w1 * w2,
    data.long = data.long,
    data.surv = data.surv,
    id = "id",
    time = "time"
  )

  # Interaction creates 3 terms: x1, x2, x1:x2
  expect_equal(ncol(result[[1]]$longitudinal$covariates), 3)
  expect_equal(ncol(result[[1]]$covariates), 3)
})

test_that("process() preserves time ordering", {
  data.long <- data.frame(
    id = c(1, 1, 1, 2, 2, 2),
    time = c(2, 0, 1, 1, 2, 0),
    v = c(6, 4, 5, 8, 9, 7),
    x = c(3, 1, 2, 5, 6, 4)
  )

  data.surv <- data.frame(
    id = c(1, 2),
    surv_time = c(3, 3),
    status = c(1, 0)
  )

  result <- process(
    formula.long = v ~ x,
    formula.surv = Surv(surv_time, status) ~ 1,
    data.long = data.long,
    data.surv = data.surv,
    id = "id",
    time = "time"
  )

  # Check ordering
  expect_equal(result[[1]]$longitudinal$times, c(0, 1, 2))
  expect_equal(unname(result[[1]]$longitudinal$measurements), c(4, 5, 6))
  expect_equal(result[[2]]$longitudinal$times, c(0, 1, 2))
  expect_equal(unname(result[[2]]$longitudinal$measurements), c(7, 8, 9))
})

test_that("process() calculates summary statistics correctly", {
  n_subjects <- 100
  n_events <- 37

  data.long <- data.frame(
    id = rep(1:n_subjects, each = 3),
    time = rep(0:2, n_subjects),
    v = rnorm(n_subjects * 3)
  )

  data.surv <- data.frame(
    id = 1:n_subjects,
    surv_time = rexp(n_subjects),
    status = c(rep(1, n_events), rep(0, n_subjects - n_events))
  )

  result <- process(
    formula.long = v ~ 1,
    formula.surv = Surv(surv_time, status) ~ 1,
    data.long = data.long,
    data.surv = data.surv,
    id = "id",
    time = "time"
  )

  expect_equal(attr(result, "n_subjects"), n_subjects)
  expect_equal(attr(result, "n_observations"), n_subjects * 3)
  expect_equal(attr(result, "event_rate"), n_events / n_subjects)
})

test_that("process() handles polynomial terms", {
  data.long <- data.frame(
    id = rep(1:10, each = 3),
    time = rep(0:2, 10),
    v = rnorm(30),
    x = runif(30)
  )

  data.surv <- data.frame(
    id = 1:10,
    surv_time = rexp(10),
    status = rbinom(10, 1, 0.5),
    z = runif(10)
  )

  result <- process(
    formula.long = v ~ poly(x, 2),
    formula.surv = Surv(surv_time, status) ~ poly(z, 3),
    data.long = data.long,
    data.surv = data.surv,
    id = "id",
    time = "time"
  )

  expect_equal(ncol(result[[1]]$longitudinal$covariates), 2)
  expect_equal(ncol(result[[1]]$covariates), 3)
})

test_that("process() handles large datasets efficiently", {
  set.seed(789)
  n_subjects <- 1000
  n_times <- 10

  data.long <- data.frame(
    id = rep(1:n_subjects, each = n_times),
    time = rep(0:(n_times - 1), n_subjects),
    v = rnorm(n_subjects * n_times),
    x1 = rnorm(n_subjects * n_times),
    x2 = rbinom(n_subjects * n_times, 1, 0.5)
  )

  data.surv <- data.frame(
    id = 1:n_subjects,
    surv_time = rexp(n_subjects, 0.1),
    status = rbinom(n_subjects, 1, 0.7),
    z = rnorm(n_subjects)
  )

  # Should complete without error
  result <- process(
    formula.long = v ~ x1 + x2,
    formula.surv = Surv(surv_time, status) ~ z,
    data.long = data.long,
    data.surv = data.surv,
    id = "id",
    time = "time"
  )

  expect_length(result, n_subjects)
  expect_equal(attr(result, "n_subjects"), n_subjects)
})

# ===========================================================================
# Additional tests for increased coverage
# ===========================================================================

test_that("process() handles subjects with no longitudinal data correctly", {
  data.long <- data.frame(
    id = c(1, 1, 3, 3), # No data for subject 2
    time = c(0, 1, 0, 1),
    v = c(1, 2, 3, 4)
  )

  data.surv <- data.frame(
    id = 1:3,
    time = c(2.5, 3.0, 2.8),
    status = c(1, 0, 1)
  )

  suppressWarnings({
    result <- process(
      v ~ 1,
      Surv(time, status) ~ 1,
      data.long,
      data.surv,
      "id",
      "time"
    )
  })

  # Check subject 2 has empty longitudinal data
  expect_equal(length(result[[2]]$longitudinal$times), 0)
  expect_equal(length(result[[2]]$longitudinal$measurements), 0)
  expect_equal(nrow(result[[2]]$longitudinal$covariates), 0)
  expect_equal(result[[2]]$longitudinal$n_obs, 0)
})

test_that("process() handles formulas with no covariates", {
  test_data <- create_minimal_test_data(
    n_subjects = 3, n_times = 3, with_covariates = FALSE
  )
  data.long <- test_data$longitudinal
  data.surv <- test_data$survival

  # Both formulas with only intercepts
  result <- process(
    v ~ 1, # No longitudinal covariates
    Surv(time, status) ~ 1, # No survival covariates
    data.long,
    data.surv,
    "id",
    "time"
  )

  # Check that covariate structures are empty
  for (i in 1:3) {
    expect_equal(ncol(result[[i]]$longitudinal$covariates), 0)
    expect_equal(ncol(result[[i]]$covariates), 0)
  }
})

test_that("process() computes event rate with all censored", {
  test_data <- create_minimal_test_data(n_subjects = 5, n_times = 2)
  data.long <- test_data$longitudinal
  data.surv <- test_data$survival
  data.surv$status <- 0 # All censored

  result <- process(
    v ~ 1,
    Surv(time, status) ~ 1,
    data.long,
    data.surv,
    "id",
    "time"
  )

  expect_equal(attr(result, "event_rate"), 0)
})

test_that("process() computes event rate with all events", {
  test_data <- create_minimal_test_data(n_subjects = 5, n_times = 2)
  data.long <- test_data$longitudinal
  data.surv <- test_data$survival
  data.surv$status <- 1 # All events

  result <- process(
    v ~ 1,
    Surv(time, status) ~ 1,
    data.long,
    data.surv,
    "id",
    "time"
  )

  expect_equal(attr(result, "event_rate"), 1)
})

test_that("process() handles edge cases in data structure", {
  # Minimum data
  test_data <- create_minimal_test_data(n_subjects = 2, n_times = 2)
  data.long.min <- test_data$longitudinal
  data.long.min$v <- c(1, 2, 3, 4) # Use specific values
  data.surv.min <- test_data$survival

  result <- process(
    v ~ 1,
    Surv(time, status) ~ 1,
    data.long.min,
    data.surv.min,
    "id",
    "time"
  )

  expect_equal(length(result), 2)
  expect_equal(attr(result, "n_subjects"), 2)
  expect_equal(attr(result, "n_observations"), 4)
})

test_that("process() handles large datasets efficiently", {
  set.seed(42)
  n_subjects <- 1000
  n_times <- 20

  data.long.large <- data.frame(
    id = rep(1:n_subjects, each = n_times),
    time = rep(seq(0, 10, length.out = n_times), n_subjects),
    v = rnorm(n_subjects * n_times),
    x1 = rnorm(n_subjects * n_times),
    x2 = rbinom(n_subjects * n_times, 1, 0.5)
  )

  data.surv.large <- data.frame(
    id = 1:n_subjects,
    time = rexp(n_subjects, 0.05) + 0.1,
    status = rbinom(n_subjects, 1, 0.6),
    w1 = rnorm(n_subjects),
    w2 = rbinom(n_subjects, 1, 0.4)
  )

  # Measure processing time
  time_taken <- system.time({
    result <- process(
      v ~ x1 + x2,
      Surv(time, status) ~ w1 + w2,
      data.long.large,
      data.surv.large,
      "id",
      "time"
    )
  })

  # Should process in reasonable time (< 5 seconds)
  expect_lt(time_taken["elapsed"], 5)
  expect_equal(length(result), n_subjects)
})
