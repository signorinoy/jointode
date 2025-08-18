# Test process function
library(survival)

test_that("process() correctly structures basic data", {
  set.seed(123)
  data.long <- data.frame(
    id = rep(1:3, each = 4),
    time = rep(0:3, 3),
    y = rnorm(12, mean = 10),
    x1 = rnorm(12),
    x2 = factor(rep(c("A", "B", "C"), each = 4))
  )

  data.surv <- data.frame(
    id = 1:3,
    obstime = c(3.5, 4.0, 3.8),
    status = c(1, 0, 1),
    z1 = c(0.5, 1.2, 0.8),
    z2 = factor(c("low", "high", "medium"))
  )

  result <- process(
    formula.long = y ~ x1 + x2,
    formula.surv = Surv(obstime, status) ~ z1 + z2,
    data.long = data.long,
    data.surv = data.surv,
    id = "id",
    time = "time"
  )

  subject1 <- result[["1"]]

  # Basic structure checks
  expect_equal(subject1$id, 1)
  expect_equal(subject1$event_time, 3.5)
  expect_equal(subject1$status, 1)

  # Longitudinal checks
  expect_equal(length(subject1$longitudinal$times), 4)
  expect_equal(subject1$longitudinal$times, 0:3)
  expect_equal(subject1$longitudinal$n_obs, 4)
  expect_true(ncol(subject1$longitudinal$covariates) >= 2)

  # Attributes
  expect_equal(attr(result, "n_subjects"), 3)
  expect_equal(attr(result, "n_observations"), 12)
  expect_equal(attr(result, "event_rate"), 2 / 3)
})

test_that("process() handles missing longitudinal data", {
  data.long <- data.frame(
    id = rep(1:2, each = 3),
    time = rep(0:2, 2),
    y = rnorm(6)
  )

  data.surv <- data.frame(
    id = 1:3,
    obstime = c(2.5, 3.0, 2.8),
    status = c(1, 0, 1)
  )

  result <- process(
    formula.long = y ~ 1,
    formula.surv = Surv(obstime, status) ~ 1,
    data.long = data.long,
    data.surv = data.surv,
    id = "id",
    time = "time"
  )

  subject3 <- result[["3"]]
  expect_equal(subject3$id, 3)
  expect_equal(subject3$event_time, 2.8)
  expect_equal(subject3$status, 1)
  expect_equal(subject3$longitudinal$n_obs, 0)
  expect_equal(length(subject3$longitudinal$times), 0)
})

test_that("process() handles formula transformations", {
  set.seed(456)
  data.long <- data.frame(
    id = rep(1:10, each = 3),
    time = rep(1:3, 10),
    y = exp(rnorm(30)),
    x = runif(30, 1, 10)
  )

  data.surv <- data.frame(
    id = 1:10,
    surv_time = rexp(10),
    status = rbinom(10, 1, 0.5),
    z = runif(10, 1, 10)
  )

  result <- process(
    formula.long = log(y) ~ log(x) + I(x^2),
    formula.surv = Surv(surv_time, status) ~ log(z),
    data.long = data.long,
    data.surv = data.surv,
    id = "id",
    time = "time"
  )

  # Check transformations applied
  expect_equal(
    unname(result[[1]]$longitudinal$measurements[1]),
    log(data.long$y[data.long$id == 1][1])
  )
  expect_equal(ncol(result[[1]]$longitudinal$covariates), 2)
})

test_that("process() handles interaction terms", {
  data.long <- data.frame(
    id = rep(1:5, each = 3),
    time = rep(0:2, 5),
    y = rnorm(15),
    x1 = rnorm(15),
    x2 = rnorm(15)
  )

  data.surv <- data.frame(
    id = 1:5,
    surv_time = rexp(5),
    status = rbinom(5, 1, 0.5),
    z1 = rnorm(5),
    z2 = rnorm(5)
  )

  result <- process(
    formula.long = y ~ x1 * x2,
    formula.surv = Surv(surv_time, status) ~ z1 * z2,
    data.long = data.long,
    data.surv = data.surv,
    id = "id",
    time = "time"
  )

  # Interaction creates 3 terms: x1, x2, x1:x2
  expect_equal(ncol(result[[1]]$longitudinal$covariates), 3)
  expect_equal(length(result[[1]]$surv_covariates), 3)
})

test_that("process() preserves time ordering", {
  data.long <- data.frame(
    id = c(1, 1, 1, 2, 2, 2),
    time = c(2, 0, 1, 1, 2, 0),
    y = c(6, 4, 5, 8, 9, 7),
    x = c(3, 1, 2, 5, 6, 4)
  )

  data.surv <- data.frame(
    id = c(1, 2),
    surv_time = c(3, 3),
    status = c(1, 0)
  )

  result <- process(
    formula.long = y ~ x,
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
    y = rnorm(n_subjects * 3)
  )

  data.surv <- data.frame(
    id = 1:n_subjects,
    surv_time = rexp(n_subjects),
    status = c(rep(1, n_events), rep(0, n_subjects - n_events))
  )

  result <- process(
    formula.long = y ~ 1,
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
    y = rnorm(30),
    x = runif(30)
  )

  data.surv <- data.frame(
    id = 1:10,
    surv_time = rexp(10),
    status = rbinom(10, 1, 0.5),
    z = runif(10)
  )

  result <- process(
    formula.long = y ~ poly(x, 2),
    formula.surv = Surv(surv_time, status) ~ poly(z, 3),
    data.long = data.long,
    data.surv = data.surv,
    id = "id",
    time = "time"
  )

  expect_equal(ncol(result[[1]]$longitudinal$covariates), 2)
  expect_equal(length(result[[1]]$surv_covariates), 3)
})

test_that("process() handles large datasets efficiently", {
  set.seed(789)
  n_subjects <- 1000
  n_times <- 10

  data.long <- data.frame(
    id = rep(1:n_subjects, each = n_times),
    time = rep(0:(n_times - 1), n_subjects),
    y = rnorm(n_subjects * n_times),
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
    formula.long = y ~ x1 + x2,
    formula.surv = Surv(surv_time, status) ~ z,
    data.long = data.long,
    data.surv = data.surv,
    id = "id",
    time = "time"
  )

  expect_length(result, n_subjects)
  expect_equal(attr(result, "n_subjects"), n_subjects)
})
