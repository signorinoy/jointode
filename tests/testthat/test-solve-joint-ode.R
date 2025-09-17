# Tests for .solve_joint_ode function

test_that(".solve_joint_ode works with basic sensitivity type", {
  # Create test data structure
  data_processed <- .process(
    longitudinal_data = sim$data$longitudinal_data,
    longitudinal_formula = observed ~ x1 + x2,
    survival_data = sim$data$survival_data,
    survival_formula = Surv(time, status) ~ w1 + w2,
    state = as.matrix(sim$data$state),
    id = "id",
    time = "time"
  )
  parameters <- sim$init

  # Test basic sensitivity type (default)
  for (i in seq_len(length(data_processed))) {
    subject_data <- data_processed[[i]]
    subject_id <- names(data_processed)[i]
    ode_solution <- .solve_joint_ode(subject_data, parameters)
    idx <- sim$data$longitudinal_data$id == subject_id
    biomarker_data <- sim$data$longitudinal_data[idx, ]
    biomarker_true <- biomarker_data$biomarker
    biomarker_pred <- ode_solution$biomarker
    names(biomarker_pred) <- NULL
    velocity_true <- biomarker_data$velocity
    velocity_pred <- ode_solution$velocity
    names(velocity_pred) <- NULL
    acceleration_true <- biomarker_data$acceleration
    acceleration_pred <- ode_solution$acceleration
    names(acceleration_pred) <- NULL

    expect_equal(biomarker_pred, biomarker_true, tolerance = 1e-2)
    expect_equal(velocity_pred, velocity_true, tolerance = 1e-3)
    expect_equal(acceleration_pred, acceleration_true, tolerance = 1e-3)
  }
})
