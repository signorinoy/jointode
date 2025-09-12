# Comprehensive validation tests

# ==============================================================================
# Formula and Data Type Validation
# ==============================================================================

test_that(".validate checks formula types", {
  # Non-formula longitudinal_formula
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = "v ~ 1",
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      id = "id",
      time = "time"
    ),
    "longitudinal_formula must be a formula"
  )

  # Non-formula survival_formula
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ 1,
      survival_formula = "Surv(time, status) ~ 1",
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      id = "id",
      time = "time"
    ),
    "survival_formula must be a formula"
  )

  # Surv() with insufficient arguments
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ 1,
      survival_formula = Surv(time) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      id = "id",
      time = "time"
    ),
    "Surv\\(\\) must have at least time and status arguments"
  )
})

test_that(".validate checks data frame types", {
  # Non-data.frame longitudinal_data
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ 1,
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = list(id = 1, time = 0, v = 1),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      id = "id",
      time = "time"
    ),
    "longitudinal_data must be a data.frame"
  )

  # Non-data.frame survival_data
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ 1,
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = list(id = 1, time = 1, status = 1),
      state = NULL,
      id = "id",
      time = "time"
    ),
    "survival_data must be a data.frame"
  )
})

test_that(".validate checks for empty data", {
  # Empty longitudinal data
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ 1,
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(
        id = numeric(0),
        time = numeric(0),
        v = numeric(0)
      ),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      id = "id",
      time = "time"
    ),
    "Longitudinal data has no rows"
  )

  # Empty survival data
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ 1,
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(
        id = numeric(0),
        time = numeric(0),
        status = numeric(0)
      ),
      state = NULL,
      id = "id",
      time = "time"
    ),
    "Survival data has no rows"
  )
})

# ==============================================================================
# Required Columns Validation
# ==============================================================================

test_that(".validate checks required columns", {
  # Missing id column in longitudinal data
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ 1,
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(time = 0, v = 1),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      id = "id",
      time = "time"
    ),
    "ID variable 'id' not found in longitudinal data"
  )

  # Missing time column in longitudinal data
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ 1,
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = 1, v = 1),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      id = "id",
      time = "time"
    ),
    "Time variable 'time' not found in longitudinal data"
  )

  # Missing id column in survival data
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ 1,
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(time = 1, status = 1),
      state = NULL,
      id = "id",
      time = "time"
    ),
    "ID variable 'id' not found in survival data"
  )
})

# ==============================================================================
# Formula Variable Validation
# ==============================================================================

test_that(".validate checks formula variables exist", {
  # Missing longitudinal response variable
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = missing_var ~ 1,
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      id = "id",
      time = "time"
    ),
    "Variables in longitudinal formula not found in data: missing_var"
  )

  # Missing longitudinal predictor variable
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ missing_pred,
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      id = "id",
      time = "time"
    ),
    "Variables in longitudinal formula not found in data: missing_pred"
  )

  # Missing survival predictor variable
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ 1,
      survival_formula = Surv(time, status) ~ missing_surv,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      id = "id",
      time = "time"
    ),
    "Variables in survival formula not found in data: missing_surv"
  )
})

# ==============================================================================
# Surv Formula Validation
# ==============================================================================

test_that(".validate checks Surv formula structure", {
  # Invalid Surv formula (no Surv on LHS)
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ 1,
      survival_formula = time ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      id = "id",
      time = "time"
    ),
    "Survival formula must have Surv\\(\\) on the left-hand side"
  )
})

# ==============================================================================
# Missing Values Validation
# ==============================================================================

test_that(".validate checks for missing values", {
  # Missing values in longitudinal data
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ 1,
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(
        id = c(1, 1),
        time = c(NA, 1),
        v = c(1, 2)
      ),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      id = "id",
      time = "time"
    ),
    "Missing values found in Time in longitudinal data"
  )

  # Missing values in survival status
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ 1,
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = NA),
      state = NULL,
      id = "id",
      time = "time"
    ),
    "Missing values in status variable"
  )
})

# ==============================================================================
# ID Consistency Validation
# ==============================================================================

test_that(".validate checks ID consistency", {
  # Duplicate IDs in survival data
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ 1,
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = rep(1, 2), time = 0:1, v = 1:2),
      survival_data = data.frame(
        id = c(1, 1),
        time = c(1, 2),
        status = c(1, 0)
      ),
      state = NULL,
      id = "id",
      time = "time"
    ),
    "Duplicate IDs found in survival data"
  )

  # Subjects in longitudinal not in survival
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ 1,
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(
        id = rep(1:3, each = 2),
        time = rep(c(0, 1), 3),
        v = 1:6
      ),
      survival_data = data.frame(id = 1:2, time = c(1, 2), status = c(1, 0)),
      state = NULL,
      id = "id",
      time = "time"
    ),
    "Subjects in longitudinal data not found in survival data: 3"
  )

  # Subjects in survival without longitudinal data (warning, not error)
  expect_warning(
    JointODE:::.validate(
      longitudinal_formula = v ~ 1,
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(
        id = rep(1:2, each = 2),
        time = rep(c(0, 1), 2),
        v = 1:4
      ),
      survival_data = data.frame(
        id = 1:4,
        time = c(1, 2, 3, 4),
        status = c(1, 0, 1, 1)
      ),
      state = NULL,
      id = "id",
      time = "time"
    ),
    "Subjects in survival data without longitudinal data: 3, 4"
  )
})

# ==============================================================================
# Time Values Validation
# ==============================================================================

test_that(".validate checks time values", {
  # Negative time in longitudinal data
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ 1,
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(
        id = c(1, 1),
        time = c(-1, 0),
        v = c(1, 2)
      ),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      id = "id",
      time = "time"
    ),
    "Negative time values found in longitudinal data"
  )

  # Non-positive time in survival data
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ 1,
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 0, status = 1),
      state = NULL,
      id = "id",
      time = "time"
    ),
    "Invalid observation times in survival data"
  )

  # Warning: measurements after observation time
  expect_warning(
    JointODE:::.validate(
      longitudinal_formula = v ~ 1,
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(
        id = rep(1:2, each = 3),
        time = rep(c(0, 1, 2), 2),
        v = 1:6
      ),
      survival_data = data.frame(
        id = 1:2,
        time = c(1.5, 2.5),
        status = c(1, 1)
      ),
      state = NULL,
      id = "id",
      time = "time"
    ),
    "1 subjects have measurements after observation time: 1"
  )

  # Warning: single longitudinal observation
  expect_warning(
    JointODE:::.validate(
      longitudinal_formula = v ~ 1,
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(
        id = 1:2,
        time = c(0, 0),
        v = c(1, 2)
      ),
      survival_data = data.frame(id = 1:2, time = c(1, 2), status = c(1, 1)),
      state = NULL,
      id = "id",
      time = "time"
    ),
    "Each subject has only one longitudinal observation"
  )
})

# ==============================================================================
# Status Values Validation
# ==============================================================================

test_that(".validate checks status values", {
  # Invalid status values
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ 1,
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 2),
      state = NULL,
      id = "id",
      time = "time"
    ),
    "Invalid status values found: 2. Must be 0 \\(censored\\) or 1 \\(event\\)"
  )
})

# ==============================================================================
# State Matrix Validation
# ==============================================================================

test_that(".validate checks state matrix", {
  # Non-matrix state
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ 1,
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = c(1, 2),
      id = "id",
      time = "time"
    ),
    "'state' must be a matrix"
  )

  # Wrong number of rows in state
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ 1,
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(
        id = rep(1:2, each = 2),
        time = rep(0:1, 2),
        v = 1:4
      ),
      survival_data = data.frame(id = 1:2, time = c(1, 2), status = c(1, 0)),
      state = matrix(c(1, 2), nrow = 1),
      id = "id",
      time = "time"
    ),
    "Invalid 'state': number of rows.*must match number of subjects"
  )

  # Wrong number of columns in state
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ 1,
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = matrix(c(1, 2, 3), nrow = 1),
      id = "id",
      time = "time"
    ),
    "Invalid 'state': must have exactly 2 columns.*got 3"
  )

  # Non-finite values in state
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ 1,
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = matrix(c(1, NA), nrow = 1),
      id = "id",
      time = "time"
    ),
    "Invalid 'state': all values must be finite"
  )
})

# ==============================================================================
# Spline Baseline Validation
# ==============================================================================

test_that(".validate checks spline_baseline parameters", {
  # Invalid parameter names
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ 1,
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      id = "id",
      time = "time",
      spline_baseline = list(invalid_param = 1)
    ),
    "Invalid parameters in spline_baseline: invalid_param"
  )

  # Invalid degree
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ 1,
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      id = "id",
      time = "time",
      spline_baseline = list(degree = 10)
    ),
    "spline_baseline\\$degree must be a single integer between 1 and 5"
  )

  # Invalid n_knots
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ 1,
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      id = "id",
      time = "time",
      spline_baseline = list(n_knots = 50)
    ),
    "spline_baseline\\$n_knots must be a single integer between 0 and 20"
  )

  # Invalid boundary_knots
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ 1,
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      id = "id",
      time = "time",
      spline_baseline = list(boundary_knots = c(0, 1, 2))
    ),
    "spline_baseline\\$boundary_knots must be NULL or.*numeric vector"
  )

  # Invalid knot_placement
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ 1,
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      id = "id",
      time = "time",
      spline_baseline = list(knot_placement = "invalid")
    ),
    "spline_baseline\\$knot_placement must be one of: quantile, equal"
  )

  # Invalid boundary_knots
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ 1,
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      id = "id",
      time = "time",
      spline_baseline = list(boundary_knots = c(2, 1))
    ),
    "boundary_knots\\[1\\] must be less than"
  )
})

# ==============================================================================
# Init Parameter Validation
# ==============================================================================

test_that(".validate checks init parameter structure", {
  # Non-list init
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ 1,
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      id = "id",
      time = "time",
      init = "not a list"
    ),
    "Invalid 'init' parameter: must be a list"
  )

  # Unknown init components
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ 1,
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      id = "id",
      time = "time",
      init = list(unknown_component = 1)
    ),
    "Invalid 'init': unknown components 'unknown_component'"
  )
})

test_that(".validate checks init$coefficients structure", {
  # Non-list coefficients
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ 1,
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      id = "id",
      time = "time",
      init = list(coefficients = "not a list")
    ),
    "Invalid 'init\\$coefficients': must be a list"
  )

  # Unknown coefficient types
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ 1,
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      id = "id",
      time = "time",
      init = list(coefficients = list(unknown_type = 1))
    ),
    "Invalid 'init\\$coefficients': unknown types 'unknown_type'"
  )
})

test_that(".validate checks init$coefficients$baseline", {
  # Non-numeric baseline
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ 1,
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      id = "id",
      time = "time",
      init = list(coefficients = list(baseline = "not numeric"))
    ),
    "Invalid 'init\\$coefficients\\$baseline': must be numeric"
  )

  # Non-finite baseline values
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ 1,
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      id = "id",
      time = "time",
      init = list(coefficients = list(baseline = c(1, NA, 3)))
    ),
    "Invalid 'init\\$coefficients\\$baseline'.*must contain finite values"
  )

  # Wrong baseline length (requires spline_baseline setup)
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ 1,
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(
        id = rep(1:2, each = 2),
        time = rep(0:1, 2),
        v = 1:4
      ),
      survival_data = data.frame(id = 1:2, time = c(2, 3), status = c(1, 0)),
      state = NULL,
      id = "id",
      time = "time",
      spline_baseline = list(degree = 3, n_knots = 2),
      init = list(coefficients = list(baseline = c(1, 2))) # Wrong length
    ),
    "Invalid 'init\\$coefficients\\$baseline'.*wrong length"
  )
})

test_that(".validate checks init$coefficients$hazard", {
  # Non-numeric hazard
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ 1,
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      id = "id",
      time = "time",
      init = list(coefficients = list(hazard = "not numeric"))
    ),
    "Invalid 'init\\$coefficients\\$hazard': must be numeric"
  )

  # Non-finite hazard values
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ 1,
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      id = "id",
      time = "time",
      init = list(coefficients = list(hazard = c(1, Inf)))
    ),
    "Invalid 'init\\$coefficients\\$hazard': must contain finite values"
  )

  # Wrong hazard length (too short)
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ x1,
      survival_formula = Surv(time, status) ~ w1,
      longitudinal_data = data.frame(
        id = c(1, 1),
        time = c(0, 1),
        v = c(1, 2),
        x1 = c(1, 2)
      ),
      survival_data = data.frame(id = 1, time = 1, status = 1, w1 = 1),
      state = NULL,
      id = "id",
      time = "time",
      init = list(coefficients = list(hazard = c(1)))
      # Should be length 3: biomarker, velocity, w1
    ),
    "Invalid 'init\\$coefficients\\$hazard'.*must have at least 2 elements"
  )

  # Wrong hazard length (exact check)
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ x1,
      survival_formula = Surv(time, status) ~ w1,
      longitudinal_data = data.frame(
        id = c(1, 1),
        time = c(0, 1),
        v = c(1, 2),
        x1 = c(1, 2)
      ),
      survival_data = data.frame(id = 1, time = 1, status = 1, w1 = 1),
      state = NULL,
      id = "id",
      time = "time",
      init = list(coefficients = list(hazard = c(1, 2, 3, 4)))
      # Should be 3, not 4
    ),
    "Invalid 'init\\$coefficients\\$hazard'.*wrong length"
  )
})

test_that(".validate checks init$coefficients$acceleration", {
  # Non-numeric acceleration
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ 1,
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      id = "id",
      time = "time",
      init = list(coefficients = list(acceleration = "not numeric"))
    ),
    "Invalid 'init\\$coefficients\\$acceleration'.*must be numeric"
  )

  # Non-finite acceleration values
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ 1,
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      id = "id",
      time = "time",
      init = list(coefficients = list(acceleration = c(1, NA)))
    ),
    "Invalid 'init\\$coefficients\\$acceleration'.*must contain finite values"
  )

  # Wrong acceleration length for autonomous model
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ x1,
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(
        id = c(1, 1),
        time = c(0, 1),
        v = c(1, 2),
        x1 = c(1, 2)
      ),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      id = "id",
      time = "time",
      autonomous = TRUE,
      init = list(coefficients = list(acceleration = c(1, 2)))
      # Should be 4 for autonomous with x1
    ),
    "Invalid 'init\\$coefficients\\$acceleration'.*wrong length.*expected 4"
  )

  # Wrong acceleration length for non-autonomous model
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ x1,
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(
        id = c(1, 1),
        time = c(0, 1),
        v = c(1, 2),
        x1 = c(1, 2)
      ),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      id = "id",
      time = "time",
      autonomous = FALSE,
      init = list(coefficients = list(acceleration = c(1, 2)))
      # Should be 5 for non-autonomous with x1
    ),
    "Invalid 'init\\$coefficients\\$acceleration'.*must have at least"
  )
})

test_that(".validate checks init$coefficients$measurement_error_sd", {
  # Non-numeric measurement_error_sd
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ 1,
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      id = "id",
      time = "time",
      init = list(coefficients = list(measurement_error_sd = "not numeric"))
    ),
    "Invalid 'init\\$coefficients\\$measurement_error_sd'.*must be"
  )

  # Non-finite measurement_error_sd
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ 1,
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      id = "id",
      time = "time",
      init = list(coefficients = list(measurement_error_sd = Inf))
    ),
    "Invalid 'init\\$coefficients\\$measurement_error_sd'.*must be finite"
  )

  # Non-positive measurement_error_sd
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ 1,
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      id = "id",
      time = "time",
      init = list(coefficients = list(measurement_error_sd = -1))
    ),
    "Invalid 'init\\$coefficients\\$measurement_error_sd'.*must be positive"
  )
})

test_that(".validate checks init$coefficients$random_effect_sd", {
  # Non-numeric random_effect_sd
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ 1,
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      id = "id",
      time = "time",
      init = list(coefficients = list(random_effect_sd = "not numeric"))
    ),
    "Invalid 'init\\$coefficients\\$random_effect_sd'.*must be"
  )

  # Non-finite random_effect_sd
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ 1,
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      id = "id",
      time = "time",
      init = list(coefficients = list(random_effect_sd = Inf))
    ),
    "Invalid 'init\\$coefficients\\$random_effect_sd': must be finite"
  )

  # Non-positive random_effect_sd
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ 1,
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      id = "id",
      time = "time",
      init = list(coefficients = list(random_effect_sd = 0))
    ),
    "Invalid 'init\\$coefficients\\$random_effect_sd': must be positive"
  )
})

test_that(".validate checks init$configurations", {
  # Non-list configurations
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ 1,
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      id = "id",
      time = "time",
      init = list(configurations = "not a list")
    ),
    "Invalid 'init\\$configurations': must be a list"
  )

  # Non-list baseline configuration
  expect_error(
    JointODE:::.validate(
      longitudinal_formula = v ~ 1,
      survival_formula = Surv(time, status) ~ 1,
      longitudinal_data = data.frame(id = c(1, 1), time = c(0, 1), v = c(1, 2)),
      survival_data = data.frame(id = 1, time = 1, status = 1),
      state = NULL,
      id = "id",
      time = "time",
      init = list(configurations = list(baseline = "not a list"))
    ),
    "Invalid 'init\\$configurations\\$baseline': must be a list"
  )
})

# ==============================================================================
# Valid Input Test
# ==============================================================================

test_that(".validate accepts valid inputs", {
  # Test passes without error or warnings
  result <- JointODE:::.validate(
    longitudinal_formula = v ~ 1,
    survival_formula = Surv(time, status) ~ 1,
    longitudinal_data = data.frame(
      id = rep(1:3, each = 2),
      time = rep(0:1, 3),
      v = 1:6
    ),
    survival_data = data.frame(
      id = 1:3,
      time = c(2, 3, 4),
      status = c(0, 1, 1)
    ),
    state = NULL,
    id = "id",
    time = "time"
  )
  expect_true(is.null(result))
})
