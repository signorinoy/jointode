# ==============================================================================
# Tests for Data Validation Functions (.validate)
# ==============================================================================
#
# Test Coverage:
#
# 1. Input Type Validation
#    - Formula objects (must be formula class)
#    - Data frames (longitudinal and survival)
#    - Character strings (id and time parameters)
#    - List objects (spline configurations)
#
# 2. Column Existence Checks
#    - ID column in both datasets
#    - Time column in longitudinal data
#    - Formula variables in respective datasets
#    - Survival time and status columns
#
# 3. Formula Structure Validation
#    - Surv() requirement for survival formula
#    - Variable availability in data
#    - Formula parsing and interpretation
#
# 4. Subject Consistency
#    - Matching subjects between datasets
#    - Orphaned longitudinal data detection
#    - Missing longitudinal data warnings
#    - Duplicate survival records
#
# 5. Temporal Consistency
#    - Non-negative longitudinal times
#    - Positive survival times
#    - Observations after event/censoring time
#
# 6. Data Integrity
#    - Status values (0 or 1 only)
#    - Missing values in critical columns
#    - Empty datasets
#    - Single observation per subject warning
#
# 7. Spline Configuration Validation
#    - Parameter names validation
#    - Degree range (1-5)
#    - Number of knots range (1-20)
#    - Knot placement strategies
#    - Boundary knot ordering
#
# 8. Edge Cases and Special Scenarios
#    - Different ID types (integer, character)
#    - Special characters in variable names
#    - Large datasets performance
#    - Completely valid data scenarios
# ==============================================================================

# ==============================================================================
# SECTION 1: INPUT TYPE VALIDATION
# ==============================================================================

test_that(".validate checks all input types correctly", {
  # Create valid test data
  long_data <- create_test_longitudinal_data(
    n_subjects = 3,
    n_times = 3,
    seed = 123
  )[, c("id", "time", "v", "x1")]

  surv_data <- create_test_survival_data(
    n_subjects = 3,
    seed = 123
  )[, c("id", "time", "status")]
  surv_data$w1 <- rnorm(3)

  long_formula <- v ~ x1
  surv_formula <- Surv(time, status) ~ w1

  # Test valid inputs pass without error
  expect_silent(
    JointODE:::.validate(
      longitudinal_formula = long_formula,
      longitudinal_data = long_data,
      survival_formula = surv_formula,
      survival_data = surv_data,
      id = "id",
      time = "time"
    )
  )

  # Test invalid longitudinal formula
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = "not a formula",
      longitudinal_data = long_data,
      survival_formula = surv_formula,
      survival_data = surv_data,
      id = "id",
      time = "time"
    ),
    "longitudinal_formula must be a formula"
  )

  # Test invalid survival formula
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = long_formula,
      longitudinal_data = long_data,
      survival_formula = "not a formula",
      survival_data = surv_data,
      id = "id",
      time = "time"
    ),
    "survival_formula must be a formula"
  )

  # Test invalid data frame
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = long_formula,
      longitudinal_data = matrix(1:12, 3, 4),
      survival_formula = surv_formula,
      survival_data = surv_data,
      id = "id",
      time = "time"
    ),
    "longitudinal_data must be a data frame"
  )

  # Test invalid id parameter (multiple values)
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = long_formula,
      longitudinal_data = long_data,
      survival_formula = surv_formula,
      survival_data = surv_data,
      id = c("id1", "id2"),
      time = "time"
    ),
    "id must be a single character string"
  )

  # Test invalid id parameter (non-character)
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = long_formula,
      longitudinal_data = long_data,
      survival_formula = surv_formula,
      survival_data = surv_data,
      id = 123,
      time = "time"
    ),
    "id must be a single character string"
  )

  # Test invalid time parameter
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = long_formula,
      longitudinal_data = long_data,
      survival_formula = surv_formula,
      survival_data = surv_data,
      id = "id",
      time = c("t1", "t2")
    ),
    "time must be a single character string"
  )
})

# ==============================================================================
# SECTION 2: COLUMN EXISTENCE VALIDATION
# ==============================================================================

test_that(".validate checks for missing columns correctly", {
  test_data <- create_minimal_test_data(
    n_subjects = 3,
    n_times = 3,
    with_covariates = FALSE
  )
  long_data <- test_data$longitudinal
  surv_data <- test_data$survival

  long_formula <- v ~ 1
  surv_formula <- Surv(time, status) ~ 1

  # Test missing id column in longitudinal data
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = long_formula,
      longitudinal_data = long_data,
      survival_formula = surv_formula,
      survival_data = surv_data,
      id = "subject_id",
      time = "time"
    ),
    "ID variable 'subject_id' not found in longitudinal data"
  )

  # Test missing time column
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = long_formula,
      longitudinal_data = long_data,
      survival_formula = surv_formula,
      survival_data = surv_data,
      id = "id",
      time = "measurement_time"
    ),
    "Time variable 'measurement_time' not found in longitudinal data"
  )

  # Test missing id in survival data
  surv_data_no_id <- surv_data
  names(surv_data_no_id)[1] <- "patient_id"
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = long_formula,
      longitudinal_data = long_data,
      survival_formula = surv_formula,
      survival_data = surv_data_no_id,
      id = "id",
      time = "time"
    ),
    "ID variable 'id' not found in survival data"
  )
})

# ==============================================================================
# SECTION 3: FORMULA VALIDATION
# ==============================================================================

test_that(".validate checks formula variables correctly", {
  test_data <- create_minimal_test_data(
    n_subjects = 3,
    n_times = 3,
    with_covariates = FALSE
  )
  long_data <- test_data$longitudinal
  surv_data <- test_data$survival

  # Test missing longitudinal formula variable
  expect_error(
    JointODE:::.validate(
      v ~ x1,
      long_data,
      Surv(time, status) ~ 1,
      surv_data,
      "id",
      "time"
    ),
    "Variables in longitudinal formula not found in data: x1"
  )

  # Test missing survival formula variable
  expect_error(
    JointODE:::.validate(
      v ~ 1,
      long_data,
      Surv(time, status) ~ w1,
      surv_data,
      "id",
      "time"
    ),
    "Variables in survival formula not found in data: w1"
  )
})

test_that(".validate checks survival formula structure", {
  test_data <- create_minimal_test_data(
    n_subjects = 3,
    n_times = 3,
    with_covariates = FALSE
  )
  long_data <- test_data$longitudinal
  surv_data <- test_data$survival

  # Test non-Surv formula
  expect_error(
    JointODE:::.validate(
      v ~ 1,
      long_data,
      time ~ status,
      surv_data,
      "id",
      "time"
    ),
    "Survival formula must have Surv\\(\\) on the left-hand side"
  )

  # Test Surv() with missing variables
  expect_error(
    JointODE:::.validate(
      v ~ 1,
      long_data,
      Surv(time2, event) ~ 1,
      surv_data,
      "id",
      "time"
    ),
    "Variables in survival formula not found in data: time2, event"
  )
})

# ==============================================================================
# SECTION 4: SUBJECT CONSISTENCY
# ==============================================================================

test_that(".validate checks subject consistency correctly", {
  # Create longitudinal data with IDs 1-4
  long_data <- data.frame(
    id = rep(1:4, each = 3),
    time = rep(0:2, 4),
    v = rnorm(12)
  )

  # Create survival data with IDs 1-3 only
  surv_data <- data.frame(
    id = 1:3,
    time = c(2.5, 3.0, 2.8),
    status = c(1, 0, 1)
  )

  # Test orphaned longitudinal data
  expect_error(
    JointODE:::.validate(
      v ~ 1,
      long_data,
      Surv(time, status) ~ 1,
      surv_data,
      "id",
      "time"
    ),
    "Subjects in longitudinal data not found in survival data: 4"
  )

  # Test warning for missing longitudinal data
  long_data_partial <- data.frame(
    id = rep(1:2, each = 3),
    time = rep(0:2, 2),
    v = rnorm(6)
  )

  expect_warning(
    JointODE:::.validate(
      v ~ 1,
      long_data_partial,
      Surv(time, status) ~ 1,
      surv_data,
      "id",
      "time"
    ),
    "Subjects in survival data without longitudinal data: 3"
  )
})

test_that(".validate handles more than 5 orphaned IDs in error message", {
  long_data_many <- data.frame(
    id = rep(1:10, each = 2),
    time = rep(0:1, 10),
    v = rnorm(20)
  )

  surv_data_few <- data.frame(
    id = 1:3,
    time = c(2, 2.5, 3),
    status = c(1, 0, 1)
  )

  expect_error(
    JointODE:::.validate(
      v ~ 1,
      long_data_many,
      Surv(time, status) ~ 1,
      surv_data_few,
      "id",
      "time"
    ),
    "Subjects in longitudinal data not found in survival data: 4, 5, 6, 7, 8"
  )
})

# ==============================================================================
# SECTION 5: TEMPORAL CONSISTENCY
# ==============================================================================

test_that(".validate checks temporal consistency correctly", {
  long_data <- create_test_longitudinal_data(
    n_subjects = 2,
    n_times = 4,
    seed = 123
  )[, c("id", "time", "v", "x1")]

  surv_data <- create_test_survival_data(
    n_subjects = 2,
    seed = 123
  )[, c("id", "time", "status")]
  # Set survival time before last observation
  surv_data$time <- c(1.5, 2.5)
  surv_data$w1 <- rnorm(2)

  # Modify longitudinal data to have observations beyond survival time
  long_data$time[long_data$id == 1 & long_data$time > 1.5] <- 2.0

  # The function issues a warning, not an error, for temporal consistency
  expect_warning(
    JointODE:::.validate(
      v ~ x1,
      long_data,
      Surv(time, status) ~ w1,
      surv_data,
      "id",
      "time"
    ),
    "subjects have measurements after observation time"
  )
})

test_that(".validate checks time values correctly", {
  # Test negative longitudinal times
  test_data <- create_minimal_test_data(n_subjects = 3, n_times = 3)
  long_data <- test_data$longitudinal
  long_data$time <- c(-1, 0, 1, 0, 1, 2, -0.5, 0.5, 1.5)
  surv_data <- test_data$survival

  expect_error(
    JointODE:::.validate(
      v ~ 1,
      long_data,
      Surv(time, status) ~ 1,
      surv_data,
      "id",
      "time"
    ),
    "Negative time values found in longitudinal data"
  )

  # Test non-positive survival times
  long_data2 <- data.frame(
    id = rep(1:3, each = 3),
    time = rep(0:2, 3),
    v = rnorm(9)
  )

  surv_data2 <- data.frame(
    id = 1:3,
    time = c(0, 3.0, 2.8),
    status = c(1, 0, 1)
  )

  expect_error(
    JointODE:::.validate(
      v ~ 1,
      long_data2,
      Surv(time, status) ~ 1,
      surv_data2,
      "id",
      "time"
    ),
    "Invalid observation times in survival data"
  )
})

# ==============================================================================
# SECTION 6: DATA INTEGRITY
# ==============================================================================

test_that(".validate checks status values correctly", {
  test_data <- create_minimal_test_data(n_subjects = 3, n_times = 3)
  long_data <- test_data$longitudinal

  # Test invalid status values
  surv_data <- test_data$survival
  surv_data$status <- c(1, 2, 0)

  expect_error(
    JointODE:::.validate(
      v ~ 1,
      long_data,
      Surv(time, status) ~ 1,
      surv_data,
      "id",
      "time"
    ),
    "Invalid status values found: 2"
  )

  # Test missing status values
  surv_data2 <- data.frame(
    id = 1:3,
    time = c(2.5, 3.0, 2.8),
    status = c(1, NA, 0)
  )

  expect_error(
    JointODE:::.validate(
      v ~ 1,
      long_data,
      Surv(time, status) ~ 1,
      surv_data2,
      "id",
      "time"
    ),
    "Missing values in status variable"
  )
})

test_that(".validate checks for NA values in critical columns", {
  # Test NA in longitudinal ID
  long_data <- data.frame(
    id = c(1, NA, 1, 2, 2, 2),
    time = c(0, 1, 2, 0, 1, 2),
    v = rnorm(6)
  )

  surv_data <- data.frame(
    id = 1:2,
    time = c(2.5, 3.0),
    status = c(1, 0)
  )

  expect_error(
    JointODE:::.validate(
      v ~ 1,
      long_data,
      Surv(time, status) ~ 1,
      surv_data,
      "id",
      "time"
    ),
    "Missing values found in ID in longitudinal data"
  )

  # Test NA in time variable
  long_data2 <- data.frame(
    id = rep(1:2, each = 3),
    time = c(0, NA, 2, 0, 1, 2),
    v = rnorm(6)
  )

  expect_error(
    JointODE:::.validate(
      v ~ 1,
      long_data2,
      Surv(time, status) ~ 1,
      surv_data,
      "id",
      "time"
    ),
    "Missing values found in Time in longitudinal data"
  )
})

test_that(".validate checks for duplicate survival records", {
  long_data <- data.frame(
    id = rep(1:3, each = 3),
    time = rep(0:2, 3),
    v = rnorm(9)
  )

  surv_data <- data.frame(
    id = c(1, 2, 2, 3),
    time = c(2.5, 3.0, 2.8, 3.2),
    status = c(1, 0, 1, 1)
  )

  expect_error(
    JointODE:::.validate(
      v ~ 1,
      long_data,
      Surv(time, status) ~ 1,
      surv_data,
      "id",
      "time"
    ),
    "Duplicate IDs found in survival data"
  )
})

# ==============================================================================
# SECTION 7: DATA DIMENSIONS
# ==============================================================================

test_that(".validate checks data dimensions correctly", {
  # Test empty longitudinal data
  long_data_empty <- data.frame(
    id = numeric(0),
    time = numeric(0),
    v = numeric(0)
  )

  surv_data <- data.frame(
    id = 1:3,
    time = c(2.5, 3.0, 2.8),
    status = c(1, 0, 1)
  )

  expect_error(
    JointODE:::.validate(
      v ~ 1,
      long_data_empty,
      Surv(time, status) ~ 1,
      surv_data,
      "id",
      "time"
    ),
    "Longitudinal data has no rows"
  )

  # Test empty survival data
  long_data <- data.frame(
    id = rep(1:3, each = 3),
    time = rep(0:2, 3),
    v = rnorm(9)
  )

  surv_data_empty <- data.frame(
    id = numeric(0),
    time = numeric(0),
    status = numeric(0)
  )

  expect_error(
    JointODE:::.validate(
      v ~ 1,
      long_data,
      Surv(time, status) ~ 1,
      surv_data_empty,
      "id",
      "time"
    ),
    "Survival data has no rows"
  )
})

test_that(".validate warns about single observations per subject", {
  # Each subject has only one observation
  test_data <- create_minimal_test_data(n_subjects = 3, n_times = 1)
  long_data <- test_data$longitudinal
  surv_data <- test_data$survival

  expect_warning(
    JointODE:::.validate(
      v ~ 1,
      long_data,
      Surv(time, status) ~ 1,
      surv_data,
      "id",
      "time"
    ),
    "Each subject has only one longitudinal observation"
  )
})

# ==============================================================================
# SECTION 8: SPLINE CONFIGURATION VALIDATION
# ==============================================================================

test_that(".validate checks spline_baseline parameters correctly", {
  test_data <- create_minimal_test_data(n_subjects = 3, n_times = 3)
  long_data <- test_data$longitudinal
  surv_data <- test_data$survival

  # Test invalid parameter names
  expect_error(
    JointODE:::.validate(
      v ~ 1,
      long_data,
      Surv(time, status) ~ 1,
      surv_data,
      "id",
      "time",
      spline_baseline = list(invalid_param = 3)
    ),
    "Invalid parameters in spline_baseline: invalid_param"
  )

  # Test invalid degree values
  expect_error(
    JointODE:::.validate(
      v ~ 1,
      long_data,
      Surv(time, status) ~ 1,
      surv_data,
      "id",
      "time",
      spline_baseline = list(degree = 0)
    ),
    "spline_baseline\\$degree must be a single integer between 1 and 5"
  )

  expect_error(
    JointODE:::.validate(
      v ~ 1,
      long_data,
      Surv(time, status) ~ 1,
      surv_data,
      "id",
      "time",
      spline_baseline = list(degree = 6)
    ),
    "spline_baseline\\$degree must be a single integer between 1 and 5"
  )

  # Test invalid n_knots values
  expect_error(
    JointODE:::.validate(
      v ~ 1,
      long_data,
      Surv(time, status) ~ 1,
      surv_data,
      "id",
      "time",
      spline_baseline = list(n_knots = -1)
    ),
    "spline_baseline\\$n_knots must be a single integer between 0 and 20"
  )

  # Test invalid knot_placement values
  expect_error(
    JointODE:::.validate(
      v ~ 1,
      long_data,
      Surv(time, status) ~ 1,
      surv_data,
      "id",
      "time",
      spline_baseline = list(knot_placement = "invalid")
    ),
    "spline_baseline\\$knot_placement must be one of: quantile, equal"
  )

  # Test invalid boundary_knots
  expect_error(
    JointODE:::.validate(
      v ~ 1,
      long_data,
      Surv(time, status) ~ 1,
      surv_data,
      "id",
      "time",
      spline_baseline = list(boundary_knots = c(2, 1))
    ),
    paste(
      "spline_baseline\\$boundary_knots\\[1\\] must be",
      "less than boundary_knots\\[2\\]"
    )
  )

  # Test valid spline_baseline parameters
  expect_silent(
    JointODE:::.validate(
      v ~ 1,
      long_data,
      Surv(time, status) ~ 1,
      surv_data,
      "id",
      "time",
      spline_baseline = list(
        degree = 3,
        n_knots = 5,
        knot_placement = "quantile",
        boundary_knots = c(0, 10)
      )
    )
  )
})

test_that(".validate checks spline_baseline as non-list", {
  test_data <- create_minimal_test_data(n_subjects = 3, n_times = 3)
  long_data <- test_data$longitudinal
  surv_data <- test_data$survival

  # Test non-list spline parameters
  expect_error(
    JointODE:::.validate(
      v ~ 1,
      long_data,
      Surv(time, status) ~ 1,
      surv_data,
      "id",
      "time",
      spline_baseline = "not a list"
    ),
    "spline_baseline must be a list"
  )
})

# ==============================================================================
# SECTION 9: EDGE CASES AND SPECIAL SCENARIOS
# ==============================================================================

test_that(".validate handles different data types correctly", {
  # Integer IDs
  long_data_int <- data.frame(
    id = as.integer(rep(1:3, each = 3)),
    time = rep(0:2, 3),
    v = rnorm(9)
  )

  surv_data_int <- data.frame(
    id = as.integer(1:3),
    time = c(2.5, 3.0, 2.8),
    status = c(1, 0, 1)
  )

  expect_silent(
    JointODE:::.validate(
      v ~ 1,
      long_data_int,
      Surv(time, status) ~ 1,
      surv_data_int,
      "id",
      "time"
    )
  )

  # Character IDs
  long_data_char <- data.frame(
    id = rep(c("A", "B", "C"), each = 3),
    time = rep(0:2, 3),
    v = rnorm(9)
  )

  surv_data_char <- data.frame(
    id = c("A", "B", "C"),
    time = c(2.5, 3.0, 2.8),
    status = c(1, 0, 1)
  )

  expect_silent(
    JointODE:::.validate(
      v ~ 1,
      long_data_char,
      Surv(time, status) ~ 1,
      surv_data_char,
      "id",
      "time"
    )
  )
})

test_that(".validate handles special characters in variable names", {
  # Data with special characters in column names
  long_data_special <- data.frame(
    `subject.id` = rep(1:3, each = 3),
    `time_point` = rep(0:2, 3),
    `outcome-1` = rnorm(9),
    check.names = FALSE
  )

  surv_data_special <- data.frame(
    `subject.id` = 1:3,
    `event.time` = c(2.5, 3.0, 2.8),
    status = c(1, 0, 1),
    check.names = FALSE
  )

  expect_silent(
    JointODE:::.validate(
      `outcome-1` ~ 1,
      long_data_special,
      Surv(`event.time`, status) ~ 1,
      surv_data_special,
      "subject.id",
      "time_point"
    )
  )
})

test_that(".validate handles large datasets efficiently", {
  skip_on_cran() # Skip performance tests on CRAN

  set.seed(42)
  n_subjects <- 1000
  n_times <- 10

  # Large longitudinal dataset
  long_data_large <- data.frame(
    id = rep(1:n_subjects, each = n_times),
    time = rep(0:(n_times - 1), n_subjects),
    v = rnorm(n_subjects * n_times),
    x1 = rnorm(n_subjects * n_times)
  )

  # Large survival dataset
  surv_data_large <- data.frame(
    id = 1:n_subjects,
    time = rexp(n_subjects, rate = 0.1) + 0.1,
    status = rbinom(n_subjects, 1, 0.7),
    w1 = rnorm(n_subjects)
  )

  # Should complete within reasonable time
  time_taken <- system.time({
    suppressWarnings(
      JointODE:::.validate(
        v ~ x1,
        long_data_large,
        Surv(time, status) ~ w1,
        surv_data_large,
        "id",
        "time"
      )
    )
  })

  # Use more generous time limit for slower systems
  expect_lt(time_taken["elapsed"], 10)
})

test_that(".validate works with completely valid data", {
  # Create completely valid data
  long_data <- create_test_longitudinal_data(
    n_subjects = 5,
    n_times = 4,
    seed = 123
  )[, c("id", "time", "v", "x1")]
  surv_data <- create_test_survival_data(n_subjects = 5, seed = 123)
  # Ensure times > max(long times)
  surv_data$time <- c(3.5, 4.0, 3.8, 4.2, 3.9)
  surv_data$w1 <- rnorm(5)

  # Should run without errors or warnings
  expect_silent(
    JointODE:::.validate(
      v ~ x1,
      long_data,
      Surv(time, status) ~ w1,
      surv_data,
      "id",
      "time"
    )
  )
})

# ==============================================================================
# SECTION 10: INIT PARAMETER VALIDATION
# ==============================================================================

test_that(".validate checks init parameter structure correctly", {
  # Create test data
  test_data <- create_minimal_test_data(n_subjects = 3, n_times = 3)
  long_data <- test_data$longitudinal
  long_data$x1 <- rnorm(nrow(long_data))
  surv_data <- test_data$survival
  surv_data$w1 <- rnorm(nrow(surv_data))

  # Test non-list init
  expect_error(
    JointODE:::.validate(
      v ~ x1,
      long_data,
      Surv(time, status) ~ w1,
      surv_data,
      "id",
      "time",
      spline_baseline = list(),
      init = "not a list"
    ),
    "Invalid 'init' parameter: must be a list"
  )

  # Test unknown top-level components
  expect_error(
    JointODE:::.validate(
      v ~ x1,
      long_data,
      Surv(time, status) ~ w1,
      surv_data,
      "id",
      "time",
      spline_baseline = list(),
      init = list(unknown_component = 1)
    ),
    paste(
      "Invalid 'init': unknown components 'unknown_component'",
      "\\(valid: coefficients, configurations\\)"
    )
  )

  # Test non-list coefficients
  expect_error(
    JointODE:::.validate(
      v ~ x1,
      long_data,
      Surv(time, status) ~ w1,
      surv_data,
      "id",
      "time",
      spline_baseline = list(),
      init = list(coefficients = "not a list")
    ),
    "Invalid 'init\\$coefficients': must be a list"
  )
})

test_that(".validate checks init coefficient types correctly", {
  test_data <- create_minimal_test_data(n_subjects = 3, n_times = 3)
  long_data <- test_data$longitudinal
  long_data$x1 <- rnorm(nrow(long_data))
  surv_data <- test_data$survival
  surv_data$w1 <- rnorm(nrow(surv_data))

  # Test unknown coefficient types
  expect_error(
    JointODE:::.validate(
      v ~ x1,
      long_data,
      Surv(time, status) ~ w1,
      surv_data,
      "id",
      "time",
      spline_baseline = list(),
      init = list(coefficients = list(unknown_coef = 1))
    ),
    "Invalid 'init\\$coefficients': unknown types 'unknown_coef'"
  )

  # Test non-numeric baseline
  expect_error(
    JointODE:::.validate(
      v ~ x1,
      long_data,
      Surv(time, status) ~ w1,
      surv_data,
      "id",
      "time",
      spline_baseline = list(),
      init = list(coefficients = list(baseline = "not numeric"))
    ),
    "Invalid 'init\\$coefficients\\$baseline': must be numeric"
  )

  # Test non-finite values
  expect_error(
    JointODE:::.validate(
      v ~ x1,
      long_data,
      Surv(time, status) ~ w1,
      surv_data,
      "id",
      "time",
      spline_baseline = list(),
      init = list(coefficients = list(baseline = c(1, 2, NA, 4, 5, 6, 7, 8)))
    ),
    "Invalid 'init\\$coefficients\\$baseline': must contain finite values"
  )
})

test_that(".validate checks init parameter lengths correctly", {
  test_data <- create_minimal_test_data(n_subjects = 3, n_times = 3)
  long_data <- test_data$longitudinal
  long_data$x1 <- rnorm(nrow(long_data))
  surv_data <- test_data$survival
  surv_data$w1 <- rnorm(nrow(surv_data))

  # Wrong baseline length (default spline has 8 basis functions)
  expect_error(
    JointODE:::.validate(
      v ~ x1,
      long_data,
      Surv(time, status) ~ w1,
      surv_data,
      "id",
      "time",
      spline_baseline = list(),
      init = list(coefficients = list(baseline = rep(0, 5)))
    ),
    paste(
      "Invalid 'init\\$coefficients\\$baseline':",
      "wrong length \\(expected 9, got 5\\)"
    )
  )

  # Wrong hazard length (should be 3 + 1 survival covariate = 4)
  expect_error(
    JointODE:::.validate(
      v ~ x1,
      long_data,
      Surv(time, status) ~ w1,
      surv_data,
      "id",
      "time",
      spline_baseline = list(),
      init = list(coefficients = list(hazard = rep(0, 2)))
    ),
    paste(
      "Invalid 'init\\$coefficients\\$hazard':",
      "must have at least 3 elements for association parameters"
    )
  )

  expect_error(
    JointODE:::.validate(
      v ~ x1,
      long_data,
      Surv(time, status) ~ w1,
      surv_data,
      "id",
      "time",
      spline_baseline = list(),
      init = list(coefficients = list(hazard = rep(0, 5)))
    ),
    paste(
      "Invalid 'init\\$coefficients\\$hazard':",
      "wrong length \\(expected 4, got 5\\)"
    )
  )

  # Wrong acceleration length (should be 3 + 1 longitudinal covariate = 4)
  expect_error(
    JointODE:::.validate(
      v ~ x1,
      long_data,
      Surv(time, status) ~ w1,
      surv_data,
      "id",
      "time",
      spline_baseline = list(),
      init = list(coefficients = list(acceleration = rep(0, 2)))
    ),
    paste(
      "Invalid 'init\\$coefficients\\$acceleration':",
      "must have at least 3 elements \\(intercept \\+ ODE parameters\\)"
    )
  )

  expect_error(
    JointODE:::.validate(
      v ~ x1,
      long_data,
      Surv(time, status) ~ w1,
      surv_data,
      "id",
      "time",
      spline_baseline = list(),
      init = list(coefficients = list(acceleration = rep(0, 6)))
    ),
    paste(
      "Invalid 'init\\$coefficients\\$acceleration':",
      "wrong length \\(expected 5, got 6\\)"
    )
  )
})

test_that(".validate checks variance parameters correctly", {
  test_data <- create_minimal_test_data(n_subjects = 3, n_times = 3)
  long_data <- test_data$longitudinal
  surv_data <- test_data$survival

  # Test non-positive measurement_error_sd
  expect_error(
    JointODE:::.validate(
      v ~ 1,
      long_data,
      Surv(time, status) ~ 1,
      surv_data,
      "id",
      "time",
      spline_baseline = list(),
      init = list(coefficients = list(measurement_error_sd = 0))
    ),
    "Invalid 'init\\$coefficients\\$measurement_error_sd': must be positive"
  )

  expect_error(
    JointODE:::.validate(
      v ~ 1,
      long_data,
      Surv(time, status) ~ 1,
      surv_data,
      "id",
      "time",
      spline_baseline = list(),
      init = list(coefficients = list(measurement_error_sd = -1))
    ),
    "Invalid 'init\\$coefficients\\$measurement_error_sd': must be positive"
  )

  # Test non-scalar measurement_error_sd
  expect_error(
    JointODE:::.validate(
      v ~ 1,
      long_data,
      Surv(time, status) ~ 1,
      surv_data,
      "id",
      "time",
      spline_baseline = list(),
      init = list(coefficients = list(measurement_error_sd = c(1, 2)))
    ),
    paste(
      "Invalid 'init\\$coefficients\\$measurement_error_sd':",
      "must be a single numeric value"
    )
  )

  # Test non-positive random_effect_sd
  expect_error(
    JointODE:::.validate(
      v ~ 1,
      long_data,
      Surv(time, status) ~ 1,
      surv_data,
      "id",
      "time",
      spline_baseline = list(),
      init = list(coefficients = list(random_effect_sd = 0))
    ),
    "Invalid 'init\\$coefficients\\$random_effect_sd': must be positive"
  )

  # Test non-finite values
  expect_error(
    JointODE:::.validate(
      v ~ 1,
      long_data,
      Surv(time, status) ~ 1,
      surv_data,
      "id",
      "time",
      spline_baseline = list(),
      init = list(coefficients = list(measurement_error_sd = Inf))
    ),
    "Invalid 'init\\$coefficients\\$measurement_error_sd': must be finite"
  )
})

test_that(".validate accepts valid init parameter", {
  test_data <- create_minimal_test_data(n_subjects = 3, n_times = 3)
  long_data <- test_data$longitudinal
  long_data$x1 <- rnorm(nrow(long_data))
  surv_data <- test_data$survival
  surv_data$w1 <- rnorm(nrow(surv_data))

  # Valid init with all components
  valid_init <- list(
    coefficients = list(
      baseline = rep(0.1, 9), # 9 spline basis functions
      hazard = rep(0.2, 4), # 3 association + 1 survival covariate
      acceleration = rep(0.3, 5), # 3 ODE params + intercept + 1 covariate
      measurement_error_sd = 0.5,
      random_effect_sd = 0.8
    ),
    configurations = list(
      baseline = list(df = 8)
    )
  )

  expect_silent(
    JointODE:::.validate(
      v ~ x1,
      long_data,
      Surv(time, status) ~ w1,
      surv_data,
      "id",
      "time",
      spline_baseline = list(),
      init = valid_init
    )
  )

  # Valid init with partial components
  partial_init <- list(
    coefficients = list(
      measurement_error_sd = 0.5,
      random_effect_sd = 0.8
    )
  )

  expect_silent(
    JointODE:::.validate(
      v ~ x1,
      long_data,
      Surv(time, status) ~ w1,
      surv_data,
      "id",
      "time",
      spline_baseline = list(),
      init = partial_init
    )
  )
})
