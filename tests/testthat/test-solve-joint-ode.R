# Tests for .solve_joint_ode function

test_that(".solve_joint_ode works with basic sensitivity type", {
  # Create test data structure
  data_processed <- .process(
    longitudinal_data = sim$data$longitudinal_data,
    longitudinal_formula = sim$formulas$longitudinal,
    survival_data = sim$data$survival_data,
    survival_formula = sim$formulas$survival,
    id = "id",
    time = "time"
  )
  parameters <- sim$parameters

  # Test basic sensitivity type (default)
  for (i in seq_len(length(data_processed))) {
    subject_data <- data_processed[[i]]
    subject_id <- names(data_processed)[i]
    ode_solution <- .solve_joint_ode(subject_data, parameters)
    idx <- sim$data$longitudinal_data$id == subject_id
    biomarker_data <- sim$data$longitudinal_data[idx, ]
    biomarker_true <- biomarker_data$biomarker
    velocity_true <- biomarker_data$velocity
    acceleration_true <- biomarker_data$acceleration
    expect_equal(ode_solution$biomarker, biomarker_true, tolerance = 1e-2)
    expect_equal(ode_solution$velocity, velocity_true, tolerance = 1e-3)
    expect_equal(ode_solution$acceleration, acceleration_true, tolerance = 1e-3)
  }
})
