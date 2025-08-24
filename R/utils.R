# Utility Functions for JointODE Package
# Internal utility functions used throughout the JointODE package.

# ===== SECTION 1: DATA VALIDATION AND PROCESSING =====

# Internal function: Validate input data for joint modeling
# Checks data integrity, formula validity, variable consistency, and temporal
# alignment between longitudinal measurements and survival outcomes
.validate <- function(
    longitudinal_formula, longitudinal_data, survival_formula, survival_data,
    id, time, spline_baseline = list(), spline_index = list()) {
  # === Step 1: Validate input types ===
  stopifnot(
    "longitudinal_formula must be a formula" =
      inherits(longitudinal_formula, "formula"),
    "longitudinal_data must be a data frame" =
      is.data.frame(longitudinal_data),
    "survival_formula must be a formula" =
      inherits(survival_formula, "formula"),
    "survival_data must be a data frame" =
      is.data.frame(survival_data),
    "id must be a single character string" =
      is.character(id) && length(id) == 1,
    "time must be a single character string" =
      is.character(time) && length(time) == 1,
    "spline_baseline must be a list" = is.list(spline_baseline),
    "spline_index must be a list" = is.list(spline_index)
  )

  # === Step 2: Validate data dimensions ===
  stopifnot(
    "Longitudinal data has no rows" = nrow(longitudinal_data) > 0,
    "Survival data has no rows" = nrow(survival_data) > 0
  )

  # === Step 3: Validate key variables exist ===
  # Check ID and time in longitudinal data
  if (!id %in% names(longitudinal_data)) {
    stop(sprintf("ID variable '%s' not found in longitudinal data", id))
  }
  if (!time %in% names(longitudinal_data)) {
    stop(sprintf("Time variable '%s' not found in longitudinal data", time))
  }

  # Check ID in survival data
  if (!id %in% names(survival_data)) {
    stop(sprintf("ID variable '%s' not found in survival data", id))
  }

  # === Step 4: Parse and validate survival formula ===
  surv_response <- survival_formula[[2]]
  if (!inherits(surv_response, "call") ||
        !identical(as.character(surv_response[[1]]), "Surv")) {
    stop("Survival formula must have Surv() on the left-hand side")
  }

  surv_call_vars <- all.vars(surv_response)
  if (length(surv_call_vars) < 2) {
    stop("Surv() must have at least time and status arguments")
  }

  time_var <- surv_call_vars[1]
  status_var <- surv_call_vars[2]

  # === Step 5: Validate formula variables exist in data ===
  # Check longitudinal formula variables
  long_vars <- all.vars(longitudinal_formula)
  missing_long_vars <- setdiff(long_vars, names(longitudinal_data))
  if (length(missing_long_vars) > 0) {
    stop(sprintf(
      "Variables in longitudinal formula not found in data: %s",
      paste(missing_long_vars, collapse = ", ")
    ))
  }

  # Check survival formula variables
  surv_vars <- all.vars(survival_formula)
  missing_surv_vars <- setdiff(surv_vars, names(survival_data))
  if (length(missing_surv_vars) > 0) {
    stop(sprintf(
      "Variables in survival formula not found in data: %s",
      paste(missing_surv_vars, collapse = ", ")
    ))
  }

  # === Step 6: Validate data integrity ===

  # Check for missing values in critical columns
  critical_cols <- list(
    "ID in longitudinal data" = longitudinal_data[[id]],
    "Time in longitudinal data" = longitudinal_data[[time]],
    "ID in survival data" = survival_data[[id]]
  )

  for (col_name in names(critical_cols)) {
    if (any(is.na(critical_cols[[col_name]]))) {
      stop(sprintf("Missing values found in %s", col_name))
    }
  }

  # Check for duplicate survival records
  if (any(duplicated(survival_data[[id]]))) {
    stop(paste(
      "Duplicate IDs found in survival data -",
      "each subject should have one record"
    ))
  }

  # === Step 7: Validate ID consistency ===
  long_ids <- unique(longitudinal_data[[id]])
  surv_ids <- unique(survival_data[[id]])

  # Check for orphaned longitudinal data
  orphaned_ids <- setdiff(long_ids, surv_ids)
  if (length(orphaned_ids) > 0) {
    stop(sprintf(
      "Subjects in longitudinal data not found in survival data: %s",
      paste(head(orphaned_ids, 5), collapse = ", ")
    ))
  }

  # Warn about missing longitudinal data
  missing_long <- setdiff(surv_ids, long_ids)
  if (length(missing_long) > 0) {
    warning(sprintf(
      "Subjects in survival data without longitudinal data: %s",
      paste(head(missing_long, 5), collapse = ", ")
    ))
  }

  # === Step 8: Validate time values ===

  # Check longitudinal times are non-negative
  if (any(longitudinal_data[[time]] < 0, na.rm = TRUE)) {
    stop("Negative time values found in longitudinal data")
  }

  # Check survival times are positive
  if (any(survival_data[[time_var]] <= 0, na.rm = TRUE)) {
    stop("Invalid observation times in survival data (must be positive)")
  }

  # Check for measurements beyond observation time
  subjects_with_late_obs <- character()
  for (subject_id in long_ids) {
    if (subject_id %in% surv_ids) {
      long_times <-
        longitudinal_data[longitudinal_data[[id]] == subject_id, time]
      surv_time <- survival_data[survival_data[[id]] == subject_id, time_var]

      if (length(surv_time) > 0 && any(long_times > surv_time + 1e-6)) {
        subjects_with_late_obs <- c(subjects_with_late_obs, subject_id)
      }
    }
  }

  if (length(subjects_with_late_obs) > 0) {
    warning(sprintf(
      "%d subjects have measurements after observation time: %s",
      length(subjects_with_late_obs),
      paste(head(subjects_with_late_obs, 3), collapse = ", ")
    ))
  }

  # === Step 9: Validate status values ===
  unique_status <- unique(survival_data[[status_var]])
  invalid_status <- setdiff(unique_status, c(0, 1, NA))

  if (length(invalid_status) > 0) {
    stop(sprintf(
      "Invalid status values found: %s. Must be 0 (censored) or 1 (event)",
      paste(invalid_status, collapse = ", ")
    ))
  }

  if (any(is.na(survival_data[[status_var]]))) {
    stop("Missing values in status variable")
  }

  # === Step 10: Check observation frequency ===
  observations_per_subject <- table(longitudinal_data[[id]])
  if (all(observations_per_subject == 1)) {
    warning(paste(
      "Each subject has only one longitudinal observation -",
      "joint modeling may not be appropriate"
    ))
  }

  # === Step 11: Validate spline parameters ===
  # Validate spline_baseline
  valid_baseline_params <- c(
    "degree", "n_knots", "knot_placement", "boundary_knots"
  )
  invalid_baseline <- setdiff(names(spline_baseline), valid_baseline_params)
  if (length(invalid_baseline) > 0) {
    stop(sprintf(
      "Invalid parameters in spline_baseline: %s. Valid parameters are: %s",
      paste(invalid_baseline, collapse = ", "),
      paste(valid_baseline_params, collapse = ", ")
    ))
  }

  if ("degree" %in% names(spline_baseline)) {
    if (!is.numeric(spline_baseline$degree) ||
          length(spline_baseline$degree) != 1 ||
          spline_baseline$degree < 1 || spline_baseline$degree > 5) {
      stop("spline_baseline$degree must be a single integer between 1 and 5")
    }
  }

  if ("n_knots" %in% names(spline_baseline)) {
    if (!is.numeric(spline_baseline$n_knots) ||
          length(spline_baseline$n_knots) != 1 ||
          spline_baseline$n_knots < 1 || spline_baseline$n_knots > 20) {
      stop("spline_baseline$n_knots must be a single integer between 1 and 20")
    }
  }

  if ("knot_placement" %in% names(spline_baseline)) {
    valid_placements <- c("quantile", "equal")
    if (!spline_baseline$knot_placement %in% valid_placements) {
      stop(sprintf(
        "spline_baseline$knot_placement must be one of: %s",
        paste(valid_placements, collapse = ", ")
      ))
    }
  }

  if ("boundary_knots" %in% names(spline_baseline)) {
    if (!is.null(spline_baseline$boundary_knots)) {
      if (
        !is.numeric(spline_baseline$boundary_knots) ||
          length(spline_baseline$boundary_knots) != 2
      ) {
        stop(paste(
          "spline_baseline$boundary_knots must be NULL or",
          "a numeric vector of length 2"
        ))
      }
      if (
        spline_baseline$boundary_knots[1] >= spline_baseline$boundary_knots[2]
      ) {
        stop(paste(
          "spline_baseline$boundary_knots[1] must be less than",
          "boundary_knots[2]"
        ))
      }
    }
  }

  # Validate spline_index
  valid_index_params <- c(
    "degree", "n_knots", "knot_placement", "boundary_knots"
  )
  invalid_index <- setdiff(names(spline_index), valid_index_params)
  if (length(invalid_index) > 0) {
    stop(sprintf(
      "Invalid parameters in spline_index: %s. Valid parameters are: %s",
      paste(invalid_index, collapse = ", "),
      paste(valid_index_params, collapse = ", ")
    ))
  }

  if ("degree" %in% names(spline_index)) {
    if (
      !is.numeric(spline_index$degree) || length(spline_index$degree) != 1 ||
        spline_index$degree < 1 || spline_index$degree > 5
    ) {
      stop("spline_index$degree must be a single integer between 1 and 5")
    }
  }

  if ("n_knots" %in% names(spline_index)) {
    if (
      !is.numeric(spline_index$n_knots) || length(spline_index$n_knots) != 1 ||
        spline_index$n_knots < 1 || spline_index$n_knots > 20
    ) {
      stop("spline_index$n_knots must be a single integer between 1 and 20")
    }
  }

  if ("knot_placement" %in% names(spline_index)) {
    valid_placements <- c("quantile", "equal")
    if (!spline_index$knot_placement %in% valid_placements) {
      stop(sprintf(
        "spline_index$knot_placement must be one of: %s",
        paste(valid_placements, collapse = ", ")
      ))
    }
  }

  if ("boundary_knots" %in% names(spline_index)) {
    if (!is.null(spline_index$boundary_knots)) {
      if (
        !is.numeric(spline_index$boundary_knots) ||
          length(spline_index$boundary_knots) != 2
      ) {
        stop(paste(
          "spline_index$boundary_knots must be NULL or",
          "a numeric vector of length 2"
        ))
      }
      if (spline_index$boundary_knots[1] >= spline_index$boundary_knots[2]) {
        stop(paste(
          "spline_index$boundary_knots[1] must be less than",
          "boundary_knots[2]"
        ))
      }
    }
  }

  invisible(NULL)
}

# Internal function: Process data for joint modeling
# Transforms longitudinal and survival data into subject-oriented structure
.process <- function(
    longitudinal_formula, longitudinal_data, survival_formula, survival_data,
    id, time) {
  # Pre-allocate with unique IDs
  unique_ids <- unique(survival_data[[id]])
  n_subjects <- length(unique_ids)
  data_process <- vector("list", n_subjects)
  names(data_process) <- as.character(unique_ids)

  # Parse survival model once
  surv_frame <- model.frame(survival_formula, data = survival_data)
  surv_response_matrix <- model.response(surv_frame)
  has_surv_covs <- length(all.vars(survival_formula[[3]])) > 0 &&
    survival_formula[[3]] != 1
  surv_design <- if (has_surv_covs) {
    model.matrix(survival_formula, surv_frame)[, -1, drop = FALSE]
  } else {
    NULL
  }


  # Build index maps for O(1) lookup
  survival_index_map <- match(unique_ids, survival_data[[id]])
  long_id_groups <- split(
    seq_len(nrow(longitudinal_data)), longitudinal_data[[id]]
  )

  for (i in seq_along(unique_ids)) {
    sid <- as.character(unique_ids[i])

    # Process longitudinal data if exists
    long_rows <- long_id_groups[[sid]]
    if (!is.null(long_rows) && length(long_rows) > 0) {
      long_subset <- longitudinal_data[long_rows, , drop = FALSE]

      # Sort once by time
      ord <- order(long_subset[[time]])
      if (!identical(ord, seq_along(ord))) {
        long_subset <- long_subset[ord, , drop = FALSE]
      }

      # Extract model components
      long_frame <- model.frame(longitudinal_formula, data = long_subset)
      long_times <- long_subset[[time]]
      long_measurements <- model.response(long_frame)
      long_covariates <- as.data.frame(
        model.matrix(longitudinal_formula, long_frame)
      )
    } else {
      long_times <- numeric(0)
      long_measurements <- numeric(0)
      long_covariates <- data.frame()
    }

    # Extract survival components (already pre-computed)
    survival_row <- survival_index_map[i]
    event_time <- surv_response_matrix[survival_row, 1]
    event_status <- surv_response_matrix[survival_row, 2]
    covariates <- if (!is.null(surv_design)) {
      surv_design[survival_row, , drop = FALSE]
    } else {
      data.frame()
    }

    # Build subject data structure
    data_process[[i]] <- list(
      id = unique_ids[i],
      time = event_time,
      status = event_status,
      covariates = covariates,
      longitudinal = list(
        times = long_times,
        measurements = long_measurements,
        covariates = long_covariates,
        n_obs = length(long_times)
      )
    )
  }

  # Compute summary statistics efficiently
  statuses <- vapply(data_process, `[[`, numeric(1), "status")
  attr(data_process, "n_subjects") <- n_subjects
  attr(data_process, "n_observations") <- nrow(longitudinal_data)
  attr(data_process, "event_rate") <- mean(statuses == 1, na.rm = TRUE)

  data_process
}

# ===== SECTION 2: NUMERICAL COMPUTATION  =====

# Configure B-Spline Parameters
# Creates B-spline configuration with optimal knot placement
.get_spline_config <- function(
    x, degree = 3, n_knots = 5, knot_placement = "quantile",
    boundary_knots = NULL) {
  # Calculate boundary knots
  if (is.null(boundary_knots)) {
    boundary_knots <- range(x, na.rm = TRUE)
  } else {
    boundary_knots <- boundary_knots
  }

  # Calculate internal knots based on placement strategy
  if (knot_placement == "quantile") {
    # Quantile-based knot placement
    probs <- seq(0, 1, length.out = n_knots + 2)[-c(1, n_knots + 2)]
    knots <- quantile(x, probs = probs, na.rm = TRUE, names = FALSE)
  } else if (knot_placement == "equal") {
    # Equally spaced knot placement
    knots <- seq(boundary_knots[1], boundary_knots[2],
      length.out = n_knots + 2
    )[-c(1, n_knots + 2)]
  } else {
    stop("knot_placement must be 'quantile' or 'equal'")
  }

  # Return configuration list
  list(
    degree = degree, knots = knots, boundary_knots = boundary_knots,
    df = length(knots) + degree + 1
  )
}

# Internal function: Compute B-spline basis using splines2
.compute_spline_basis <- function(x, config) {
  splines2::bSpline(
    x = x,
    knots = config$knots,
    degree = config$degree,
    Boundary.knots = config$boundary_knots,
    intercept = TRUE,
    warn.outside = FALSE
  )
}

# Internal function: Compute B-spline basis derivative using splines2
.compute_spline_basis_deriv <- function(x, config) {
  # splines2::dbs computes derivatives of B-splines efficiently
  splines2::dbs(
    x = x,
    knots = config$knots,
    degree = config$degree,
    Boundary.knots = config$boundary_knots,
    derivs = 1,
    intercept = TRUE,
    warn.outside = FALSE
  )
}

# ===== SECTION 3: ODE SOLVING HELPERS =====

# Internal function: Get longitudinal covariates at time point (optimized)
.get_longitudinal_covariates <- function(data, time = NULL, row_index = NULL) {
  # Fast return for empty case
  long_cov <- data$longitudinal$covariates
  if (is.null(long_cov) || nrow(long_cov) == 0) {
    return(numeric(0))
  }

  # If time is provided, find nearest time point
  if (is.null(row_index)) {
    if (!is.null(time) && !is.null(data$longitudinal$times)) {
      row_index <- which.min(abs(data$longitudinal$times - time))
    } else {
      row_index <- nrow(long_cov)
    }
  }

  as.numeric(long_cov[row_index, , drop = TRUE])
}

# Internal function: Compute acceleration function
.compute_acceleration <- function(
    biomarker, velocity, time, data, coefficients, config) {
  long_cov <- .get_longitudinal_covariates(data, time)
  z <- as.vector(c(biomarker, velocity, long_cov, time))
  index_value <- sum(coefficients$index_beta * z)

  # Compute basis and return
  basis_g <- .compute_spline_basis(index_value, config$index)
  sum(basis_g * coefficients$index_g)
}

.compute_acceleration_deriv <- function(
    biomarker, velocity, time, data, coefficients, config, type = "beta",
    dbiomarker = NULL, dvelocity = NULL) {
  # Get longitudinal covariates
  long_cov <- .get_longitudinal_covariates(data, time)

  # Construct Z vector: [m(t), m_dot(t), X(t), t]
  z_vec <- c(biomarker, velocity, long_cov, time)

  # Get parameters
  beta <- coefficients$index_beta
  theta <- coefficients$index_g

  # Compute index value u = beta' * Z
  index_value <- sum(beta * z_vec)

  # Compute B-spline basis and its derivative
  basis_g <- .compute_spline_basis(index_value, config$index)
  basis_g_deriv <- .compute_spline_basis_deriv(index_value, config$index)

  if (type == "beta") {
    # Compute d(acceleration)/dbeta
    # d(ddot{m})/dbeta = theta' * B'_g(u) * du/dbeta
    # where du/dbeta = Z + beta' * dZ/dbeta

    n_beta <- length(beta)

    if (!is.null(dbiomarker) && !is.null(dvelocity)) {
      # Full chain rule with sensitivity tracking
      dz_dbeta <- matrix(0, length(z_vec), n_beta)
      dz_dbeta[1, ] <- dbiomarker
      dz_dbeta[2, ] <- dvelocity

      # Total derivative: du/dbeta_j = z_j + beta' * dZ/dbeta_j
      # This gives a vector of length n_beta
      du_dbeta <- numeric(n_beta)
      for (j in seq_len(n_beta)) {
        if (j <= length(z_vec)) {
          du_dbeta[j] <- z_vec[j] + sum(beta * dz_dbeta[, j])
        } else {
          du_dbeta[j] <- sum(beta * dz_dbeta[, j])
        }
      }
    } else {
      # If sensitivities not provided, use just Z
      # (This would be an approximation)
      # Each beta_j affects u through z_j
      du_dbeta <- numeric(n_beta)
      for (j in seq_len(n_beta)) {
        if (j <= length(z_vec)) {
          du_dbeta[j] <- z_vec[j]
        }
      }
    }

    # Compute d(acceleration)/dbeta = theta' * B'_g(u) * du/dbeta
    dacceleration_dbeta <- sum(theta * basis_g_deriv) * du_dbeta

    return(dacceleration_dbeta)

  } else if (type == "theta") {
    # Compute d(acceleration)/dtheta
    # d(ddot{m})/dtheta = B_g(u) + theta' * B'_g(u) * du/dtheta
    # where du/dtheta = beta' * dZ/dtheta

    du_dtheta <- beta[1] * dbiomarker + beta[2] * dvelocity

    # Total derivative:
    # d(acceleration)/dtheta = B_g(u) + theta' * B'_g(u) * du/dtheta
    dacceleration_dtheta <- as.vector(basis_g) +
      sum(theta * basis_g_deriv) * du_dtheta


    return(dacceleration_dtheta)

  } else {
    stop("Invalid type. Must be one of: 'beta', 'theta'")
  }
}

# Internal function: Compute log-hazard function
.compute_log_hazard <- function(
    time, biomarker, velocity, acceleration, data, coefficients, config) {
  # Baseline hazard
  basis_lambda <- .compute_spline_basis(time, config$baseline)
  log_baseline <- sum(basis_lambda * coefficients$baseline)

  # Biomarker effects
  biomarker_vec <- c(biomarker, velocity, acceleration)
  log_biomarker <- sum(biomarker_vec * coefficients$hazard[1:3])

  # Covariate effects
  n_hazard <- length(coefficients$hazard)
  log_covariate <- if (n_hazard > 3 && !is.null(data$covariates)) {
    w <- data$covariates
    if (nrow(w) > 0) {
      sum(coefficients$hazard[4:n_hazard] * w)
    } else {
      0
    }
  } else {
    0
  }

  log_baseline + log_biomarker + log_covariate
}

# Internal function: Prepare initial conditions for ODE
.prepare_initial_conditions <- function(
    biomarker_initial, coefficients, sensitivity_type) {
  if (length(biomarker_initial) != 2) {
    stop("biomarker_initial must be a numeric vector of length 2")
  }

  switch(sensitivity_type,
    basic = c(0, biomarker_initial),
    eta_alpha = {
      n_basis <- coefficients$config$baseline$df
      c(0, biomarker_initial, rep(0, n_basis + 3))
    },
    beta = {
      n_index_beta <- length(coefficients$coef$index_beta)
      c(0, biomarker_initial, rep(0, 3 * n_index_beta))
    },
    theta = {
      n_index_basis <- length(coefficients$coef$index_g)
      c(0, biomarker_initial, rep(0, 3 * n_index_basis))
    }
  )
}

# Internal function: Extract results from ODE solution
.extract_ode_results <- function(solution, data, parameters, sensitivity_type) {
  # Get final state
  final_state <- solution[nrow(solution), -1]
  event_time <- data$time

  # Compute final values
  biomarker_final <- final_state[2]
  velocity_final <- final_state[3]

  acceleration_final <- .compute_acceleration(
    biomarker_final, velocity_final, event_time,
    data, parameters$coef, parameters$config
  )

  log_hazard_final <- .compute_log_hazard(
    event_time, biomarker_final, velocity_final,
    acceleration_final, data, parameters$coef, parameters$config
  )

  # Extract biomarker trajectory
  biomarker_values <- if (length(data$longitudinal$times) > 0) {
    obs_indices <- match(data$longitudinal$times, solution[, 1])
    solution[obs_indices, 3]
  } else {
    NULL
  }
  velocity_values <- if (length(data$longitudinal$times) > 0) {
    obs_indices <- match(data$longitudinal$times, solution[, 1])
    solution[obs_indices, 4]
  } else {
    NULL
  }
  acceleration_values <- if (length(data$longitudinal$times) > 0) {
    obs_indices <- match(data$longitudinal$times, solution[, 1])
    sapply(seq_along(obs_indices), function(i) {
      idx <- obs_indices[i]
      .compute_acceleration(
        solution[idx, 3], solution[idx, 4], solution[idx, 1],
        data, parameters$coef, parameters$config
      )
    })
  } else {
    NULL
  }

  # Build base result
  result <- list(
    cum_hazard = final_state[1],
    log_hazard = log_hazard_final,
    biomarker = biomarker_values,
    velocity = velocity_values,
    acceleration = acceleration_values
  )

  # Add sensitivity-specific outputs
  if (sensitivity_type == "eta_alpha") {
    n_basis <- parameters$config$baseline$df
    result <- list(
      log_hazard = log_hazard_final,
      cum_hazard = final_state[1],
      biomarker = final_state[2],
      velocity = final_state[3],
      acceleration = acceleration_final,
      d_cum_hazard_d_eta = final_state[4:(4 + n_basis - 1)],
      d_cum_hazard_d_alpha = final_state[(4 + n_basis):(4 + n_basis + 2)]
    )
  } else if (sensitivity_type == "beta") {
    n_beta <- length(parameters$coef$index_beta)
    dbiomarker_dbeta_values <- if (length(data$longitudinal$times) > 0) {
      obs_indices <- match(data$longitudinal$times, solution[, 1])
      solution[obs_indices, 5:(4 + n_beta), drop = FALSE]
    } else {
      NULL
    }
    dvelocity_dbeta_final <- final_state[(4 + n_beta):(3 + 2 * n_beta)]
    dacceleration_dbeta_final <- .compute_acceleration_deriv(
      biomarker_final, velocity_final, event_time,
      data, parameters$coef, parameters$config, type = "beta",
      dbiomarker = final_state[4:(3 + n_beta)],
      dvelocity = dvelocity_dbeta_final
    )
    dcum_hazard_dbeta_final <- final_state[(4 + 2 * n_beta):(3 + 3 * n_beta)]
    dbiomarker_dbeta_at_event <- final_state[4:(3 + n_beta)]
    result <- list(
      log_hazard = log_hazard_final,
      cum_hazard = final_state[1],
      biomarker = biomarker_values,
      dbiomarker_dbeta = dbiomarker_dbeta_values,
      dbiomarker_dbeta_at_event = dbiomarker_dbeta_at_event,
      dvelocity_dbeta = dvelocity_dbeta_final,
      dacceleration_dbeta = dacceleration_dbeta_final,
      dcum_hazard_dbeta = dcum_hazard_dbeta_final
    )
  } else if (sensitivity_type == "theta") {
    n_theta <- length(parameters$coef$index_g)
    return(n_theta)
  }

  result
}

# Internal function: Solve joint ODE system
# Solves the coupled ODE system for cumulative hazard and biomarker trajectory
# Supports multiple configurations: basic, eta_alpha, beta, theta
.solve_joint_ode <- function(data, parameters, sensitivity_type = "basic") {
  # Validate sensitivity_type
  valid_types <- c("basic", "eta_alpha", "beta", "theta")
  if (!sensitivity_type %in% valid_types) {
    stop(
      "Invalid sensitivity_type. Must be one of: ",
      paste(valid_types, collapse = ", ")
    )
  }

  # Define ODE derivatives function based on sensitivity type
  ode_derivatives <- function(time, state, parameters) {
    biomarker <- state[2]
    velocity <- state[3]

    # Compute acceleration
    acceleration <- .compute_acceleration(
      biomarker, velocity, time,
      parameters$data, parameters$coef, parameters$config
    )

    # Compute hazard (without random effect b)
    log_hazard <- .compute_log_hazard(
      time, biomarker, velocity, acceleration,
      parameters$data, parameters$coef, parameters$config
    )
    hazard <- exp(log_hazard)

    # Basic derivatives: [dΛ/dt, dm/dt, dṁ/dt]
    basic_derivs <- c(hazard, velocity, acceleration)

    if (parameters$sensitivity_type == "basic") {
      list(basic_derivs)
    } else if (parameters$sensitivity_type == "eta_alpha") {
      basis_lambda <- .compute_spline_basis(time, parameters$config$baseline)
      m_vec <- c(biomarker, velocity, acceleration)
      list(c(basic_derivs, as.numeric(basis_lambda) * hazard, m_vec * hazard))
    } else if (parameters$sensitivity_type == "beta") {
      n_beta <- length(parameters$coef$index_beta)
      dbiomarker_dbeta <- state[4:(3 + n_beta)]
      dvelocity_dbeta <- state[(4 + n_beta):(3 + 2 * n_beta)]
      dacceleration_dbeta <- .compute_acceleration_deriv(
        biomarker, velocity, time, data, parameters$coef, parameters$config,
        type = "beta",
        dbiomarker = dbiomarker_dbeta,
        dvelocity = dvelocity_dbeta
      )
      # Compute dLambda/dbeta derivative
      alpha <- parameters$coef$hazard[1:3]
      dm_vec_dbeta <- rbind(
        dbiomarker_dbeta,
        dvelocity_dbeta,
        dacceleration_dbeta
      )
      d_lambda_dt_dbeta <- as.vector(
        crossprod(alpha, dm_vec_dbeta)
      ) * hazard

      # Return derivatives
      list(c(
        basic_derivs, dvelocity_dbeta, dacceleration_dbeta, d_lambda_dt_dbeta
      ))
    } else if (sensitivity_type == "theta") {
      n_theta <- length(parameters$coef$index_g)
      dbiomarker_dtheta <- state[4:(3 + n_theta)]
      dvelocity_dtheta <- state[(4 + n_theta):(3 + 2 * n_theta)]
      dacceleration_dtheta <- .compute_acceleration_deriv(
        biomarker, velocity, time, data, parameters$coef, parameters$config,
        type = "theta",
        dbiomarker = dbiomarker_dtheta,
        dvelocity = dvelocity_dtheta
      )

      # Compute dLambda/dtheta derivative
      alpha <- parameters$coef$hazard[1:3]
      dm_vec_dtheta <- rbind(
        dbiomarker_dtheta,
        dvelocity_dtheta,
        dacceleration_dtheta
      )
      d_lambda_dt_dtheta <- as.vector(
        crossprod(alpha, dm_vec_dtheta)
      ) * hazard

      # Return derivatives
      list(c(
        basic_derivs, dvelocity_dtheta, dacceleration_dtheta, d_lambda_dt_dtheta
      ))
    } else {
      stop("Unknown sensitivity type: ", sensitivity_type)
    }
  }

  # Prepare initial conditions based on sensitivity type
  initial_extended <- .prepare_initial_conditions(
    c(0, 0), parameters, sensitivity_type
  )

  # Solve ODE System
  ode_parameters <- list(
    data = data,
    coef = parameters$coef,
    config = parameters$config,
    sensitivity_type = sensitivity_type
  )
  event_time <- data$time
  times <- sort(unique(c(0, data$longitudinal$times, event_time)))

  n_timepoints <- length(times)

  # Adaptive tolerances - looser for larger problems
  rtol <- if (n_timepoints > 50) 1e-2 else 1e-3
  atol <- if (n_timepoints > 50) 1e-3 else 1e-4

  # Choose optimal method

  solution <- deSolve::ode(
    y = initial_extended,
    times = times,
    func = ode_derivatives,
    parms = ode_parameters,
    method = "ode45",
    rtol = rtol,
    atol = atol
  )

  # Extract and return results
  .extract_ode_results(solution, data, parameters, sensitivity_type)
}

# ===== SECTION 4: STATISTICAL DISTRIBUTIONS =====

# Internal function: Compute posterior moments using AGHQ
# Uses adaptive Gauss-Hermite quadrature for random effect posterior
.compute_posterior_aghq <- function(
    ode_solution, data, b_hat_init, measurement_error_sd, random_effect_sd,
    k = 7) {
  # Extract components
  cum_hazard_0 <- ode_solution$cum_hazard
  biomarker_0 <- ode_solution$biomarker

  status <- data$status
  measurements <- data$longitudinal$measurements
  n_obs <- data$longitudinal$n_obs

  # Compute residual sum
  s_i <- sum(measurements - biomarker_0)
  inv_measurement_error_sd2 <- 1 / measurement_error_sd^2
  inv_random_effect_sd2 <- 1 / random_effect_sd^2

  # Define log-posterior function (up to normalizing constant)
  logpost <- function(b) {
    -b^2 * (n_obs * inv_measurement_error_sd2 + inv_random_effect_sd2) +
      b * (s_i * inv_measurement_error_sd2 + status) -
      exp(b) * cum_hazard_0
  }
  logpost_grad <- function(b) {
    2 * b * (n_obs * inv_measurement_error_sd2 + inv_random_effect_sd2) +
      (s_i * inv_measurement_error_sd2 + status) -
      exp(b) * cum_hazard_0
  }
  logpost_hess <- function(b) {
    -2 * (n_obs * inv_measurement_error_sd2 + inv_random_effect_sd2) -
      exp(b) * cum_hazard_0
  }

  # Build adaptive quadrature
  fit_aghq <- aghq::aghq(
    ff = list(
      fn = logpost,
      gr = logpost_grad,
      he = logpost_hess
    ),
    k = k,
    startingvalue = b_hat_init,
  )

  # Compute posterior moments
  b_hat <- aghq::compute_moment(
    fit_aghq$normalized_posterior,
    ff = function(x) x
  )
  b2_hat <- aghq::compute_moment(
    fit_aghq$normalized_posterior,
    ff = function(x) x^2
  )
  v_hat <- b2_hat - b_hat^2
  # Compute E[exp(b)] for survival component
  exp_b <- aghq::compute_moment(
    fit_aghq$normalized_posterior,
    ff = function(x) exp(x)
  )

  list(
    b_hat = b_hat,
    v_hat = v_hat,
    exp_b = exp_b
  )
}

# Internal function: Compute Q function for eta, alpha, phi only
# Fixes beta and theta parameters while optimizing hazard-related parameters
.compute_q_eta_alpha_phi <- function(
    params, data_list, posteriors, config, fixed_params) {
  # Get dimensions
  n_baseline <- config$baseline$df
  n_surv_covariates <- ncol(data_list[[1]]$covariates)

  # Parse parameter vector
  idx <- 1
  eta <- params[idx:(idx + n_baseline - 1)] # Baseline hazard coefficients
  idx <- idx + n_baseline
  alpha <- params[idx:(idx + 2)] # Biomarker hazard coefficients
  idx <- idx + 3
  phi <- params[idx:(idx + n_surv_covariates - 1)]

  # Build parameters with fixed beta and theta
  parameters <- list(
    coef = list(
      baseline = eta,
      hazard = c(alpha, phi),
      index_g = fixed_params$index_g, # Fixed g(u) coefficients
      index_beta = fixed_params$index_beta # Fixed index coefficients
    ),
    config = config
  )

  # Pre-extract all data for vectorized operations
  n_subjects <- length(data_list)
  subject_ids <- names(data_list)

  # Batch extract all data
  times <- vapply(data_list, `[[`, numeric(1), "time")
  status_vec <- vapply(data_list, `[[`, numeric(1), "status")
  exp_b_vec <- vapply(
    subject_ids, function(id) posteriors[[id]]$exp_b, numeric(1)
  )

  # Extract covariates as matrix
  covariates_mat <- do.call(rbind, lapply(data_list, `[[`, "covariates"))

  # Compute basis once for all subjects
  basis_lambda <- .compute_spline_basis(times, config$baseline)

  # Pre-allocate storage for ODE results
  log_hazard_vec <- numeric(n_subjects)
  cum_hazard_vec <- numeric(n_subjects)
  biomarkers_mat <- matrix(0, n_subjects, 3)
  d_cum_hazard_d_eta_mat <- matrix(0, n_subjects, n_baseline)
  d_cum_hazard_d_alpha_mat <- matrix(0, n_subjects, 3)

  # Compute all ODE solutions with optimized batch processing
  # Group subjects by number of time points for batch optimization
  n_timepoints_vec <- vapply(
    data_list, function(d) length(d$longitudinal$times), integer(1)
  )
  subject_order <- order(n_timepoints_vec)

  for (idx in seq_len(n_subjects)) {
    i <- subject_order[idx]
    ode_sol <- .solve_joint_ode(
      data_list[[i]], parameters,
      sensitivity_type = "eta_alpha"
    )
    log_hazard_vec[i] <- ode_sol$log_hazard
    cum_hazard_vec[i] <- ode_sol$cum_hazard
    biomarkers_mat[i, ] <- c(
      ode_sol$biomarker, ode_sol$velocity, ode_sol$acceleration
    )
    d_cum_hazard_d_eta_mat[i, ] <- ode_sol$d_cum_hazard_d_eta
    d_cum_hazard_d_alpha_mat[i, ] <- ode_sol$d_cum_hazard_d_alpha
  }

  q <- sum(status_vec * log_hazard_vec) - sum(exp_b_vec * cum_hazard_vec)

  # Gradient w.r.t. eta (baseline hazard)
  grad_eta <- as.vector(crossprod(status_vec, basis_lambda)) -
    as.vector(crossprod(exp_b_vec, d_cum_hazard_d_eta_mat))

  # Gradient w.r.t. alpha (biomarker effects)
  grad_alpha <- as.vector(crossprod(status_vec, biomarkers_mat)) -
    as.vector(crossprod(exp_b_vec, d_cum_hazard_d_alpha_mat))

  # Gradient w.r.t. phi (covariate effects)
  grad_phi <- as.vector(crossprod(status_vec, covariates_mat)) -
    as.vector(crossprod(exp_b_vec * cum_hazard_vec, covariates_mat))

  # Combine gradient components
  q_grad <- c(grad_eta, grad_alpha, grad_phi)


  # Return negative Q for minimization
  res <- -q
  attr(res, "gradient") <- -q_grad
  res
}

.compute_q_eta_alpha_phi_grad <- function(
    params, data_list, posteriors, config, fixed_params) {
  attr(
    .compute_q_eta_alpha_phi(
      params, data_list, posteriors, config, fixed_params
    ),
    "gradient"
  )
}

.compute_q_beta <- function(
    params, data_list, posteriors, config, fixed_params) {
  beta <- params

  # Get dimensions
  n_beta <- length(beta)
  n_subjects <- length(data_list)

  # Extract fixed parameters
  measurement_error_sd <- fixed_params$measurement_error_sd
  inv_sigma_e2 <- 1 / measurement_error_sd^2
  alpha <- fixed_params$hazard[1:3]
  theta <- fixed_params$index_g

  # Parameters for ODE solving
  parameters <- list(
    coef = list(
      baseline = fixed_params$baseline,
      hazard = fixed_params$hazard,
      index_g = theta,
      index_beta = beta
    ),
    config = config
  )

  # Extract posterior expectations
  subject_ids <- names(data_list)
  b_hat_vec <- vapply(
    subject_ids, function(id) posteriors[[id]]$b_hat, numeric(1)
  )
  exp_b_vec <- vapply(
    subject_ids, function(id) posteriors[[id]]$exp_b, numeric(1)
  )
  status_vec <- vapply(data_list, `[[`, numeric(1), "status")

  # Pre-allocate storage for ODE results
  log_hazard_vec <- numeric(n_subjects)
  cum_hazard_vec <- numeric(n_subjects)
  dcum_hazard_dbeta <- matrix(0, n_subjects, n_beta)

  residuals_list <- vector("list", n_subjects)
  dbiomarker_dbeta_list <- vector("list", n_subjects)
  dbiomarker_dbeta_final <- matrix(0, n_subjects, n_beta)
  dvelocity_dbeta_final <- matrix(0, n_subjects, n_beta)
  dacceleration_dbeta_final <- matrix(0, n_subjects, n_beta)

  # Sort subjects by number of timepoints for efficient ODE solving
  n_timepoints_vec <- vapply(
    data_list, function(d) length(d$longitudinal$times), integer(1)
  )
  subject_order <- order(n_timepoints_vec)

  # Batch solve ODEs for all subjects
  for (idx in seq_len(n_subjects)) {
    i <- subject_order[idx]
    subject_data <- data_list[[i]]

    # Solve ODE with beta sensitivities
    ode_sol <- .solve_joint_ode(
      subject_data, parameters, sensitivity_type = "beta"
    )

    # Store ODE results
    log_hazard_vec[i] <- ode_sol$log_hazard
    cum_hazard_vec[i] <- ode_sol$cum_hazard

    # Compute and store residuals
    measurements <- subject_data$longitudinal$measurements
    residuals_list[[i]] <- measurements - ode_sol$biomarker - b_hat_vec[i]

    # Store derivatives at observation times
    dbiomarker_dbeta_list[[i]] <- ode_sol$dbiomarker_dbeta

    # Extract derivatives at final time for survival contribution
    # We need derivatives at event time from the ODE final state
    dbiomarker_dbeta_final[i, ] <- ode_sol$dbiomarker_dbeta_at_event
    dvelocity_dbeta_final[i, ] <- ode_sol$dvelocity_dbeta
    dacceleration_dbeta_final[i, ] <- ode_sol$dacceleration_dbeta
    dcum_hazard_dbeta[i, ] <- ode_sol$dcum_hazard_dbeta
  }

  # Compute Q value components
  q_long <- 0
  grad_long <- numeric(n_beta)
  for (idx in seq_len(n_subjects)) {
    i <- subject_order[idx]
    residuals <- residuals_list[[i]]
    n_obs_i <- length(residuals)
    if (n_obs_i > 0) {
      q_long <- q_long - sum(residuals^2) * inv_sigma_e2 / 2
      grad_long <- grad_long + as.vector(
        crossprod(residuals, dbiomarker_dbeta_list[[i]]) * inv_sigma_e2
      )
    }
  }
  q_surv <- sum(status_vec * log_hazard_vec - exp_b_vec * cum_hazard_vec)
  grad_surv <- as.vector(
    crossprod(status_vec, dbiomarker_dbeta_final) * alpha[1] +
      crossprod(status_vec, dvelocity_dbeta_final) * alpha[2] +
      crossprod(status_vec, dacceleration_dbeta_final) * alpha[3]
  ) - as.vector(crossprod(exp_b_vec, dcum_hazard_dbeta))

  # Combine Q components
  q_total <- -(q_long + q_surv)
  attr(q_total, "gradient") <- -(grad_long + grad_surv)
  q_total
}

# Helper function: Compute gradient of Q w.r.t. beta
.compute_q_beta_grad <- function(
    params, data_list, posteriors, config, fixed_params) {
  attr(
    .compute_q_beta(
      params, data_list, posteriors, config, fixed_params
    ),
    "gradient"
  )
}

# Internal function: Compute Q function for M-step
# Computes expected complete-data log-likelihood for EM algorithm
.compute_q_function <- function(
    theta, data_list, posteriors, config,
    measurement_error_sd, random_effect_sd) {
  # Get dimensions for each parameter group
  n_baseline <- config$baseline$df
  n_index_g <- config$index$df

  # Find first subject with longitudinal data to get covariate dimensions
  subject_with_long <- NULL
  for (subject in data_list) {
    if (subject$longitudinal$n_obs > 0) {
      subject_with_long <- subject
      break
    }
  }

  # Calculate dimensions
  n_long_covariates <- if (!is.null(subject_with_long)) {
    ncol(subject_with_long$longitudinal$covariates)
  } else {
    0
  }
  n_surv_covariates <- ncol(data_list[[1]]$covariates)
  n_index_beta <- n_long_covariates + 3
  n_hazard <- n_surv_covariates + 3

  # Parse parameter vector sequentially
  idx <- 1
  baseline_coef <- theta[idx:(idx + n_baseline - 1)]
  idx <- idx + n_baseline
  hazard_coef <- theta[idx:(idx + n_hazard - 1)]
  idx <- idx + n_hazard
  index_g_coef <- theta[idx:(idx + n_index_g - 1)]
  idx <- idx + n_index_g
  index_beta_coef <- theta[idx:(idx + n_index_beta - 1)]

  parameters <- list(
    coef = list(
      baseline = baseline_coef,
      hazard = hazard_coef,
      index_g = index_g_coef,
      index_beta = index_beta_coef
    ),
    config = config
  )

  # Pre-compute constants for efficiency
  inv_2sigma2 <- 1 / (2 * measurement_error_sd^2)
  log_2pi_sigma2 <- log(2 * pi * measurement_error_sd^2)
  half_log_2pi_sigma2 <- 0.5 * log_2pi_sigma2

  # Initialize Q components
  q_long <- 0
  q_surv <- 0
  n_total_obs <- 0

  # Vectorized computation using lapply
  subject_ids <- names(data_list)

  # Pre-extract all posterior values for efficiency
  b_hat_all <- vapply(
    subject_ids, function(id) posteriors[[id]]$b_hat, numeric(1)
  )
  v_hat_all <- vapply(
    subject_ids, function(id) posteriors[[id]]$v_hat, numeric(1)
  )
  exp_b_all <- vapply(
    subject_ids, function(id) posteriors[[id]]$exp_b, numeric(1)
  )

  # Compute all ODE solutions
  ode_solutions <- lapply(data_list, function(subject_data) {
    .solve_joint_ode(subject_data, parameters)
  })

  # Process each subject and compute Q components
  subject_results <- mapply(function(subject_data, ode_solution, idx) {
    b_hat <- b_hat_all[idx]
    v_hat <- v_hat_all[idx]
    exp_b <- exp_b_all[idx]

    # Longitudinal component
    q_long_i <- 0
    n_obs_i <- 0
    if (subject_data$longitudinal$n_obs > 0) {
      residuals <- subject_data$longitudinal$measurements -
        ode_solution$biomarker - b_hat
      q_long_i <- -sum(residuals^2 + v_hat) * inv_2sigma2 -
        subject_data$longitudinal$n_obs * half_log_2pi_sigma2
      n_obs_i <- subject_data$longitudinal$n_obs
    }

    # Survival component
    q_surv_i <- -exp_b * ode_solution$cum_hazard
    if (subject_data$status == 1) {
      q_surv_i <- q_surv_i + ode_solution$log_hazard + b_hat
    }

    list(q_long = q_long_i, q_surv = q_surv_i, n_obs = n_obs_i)
  }, data_list, ode_solutions, seq_along(data_list), SIMPLIFY = FALSE)

  # Sum all components
  q_long <- sum(vapply(subject_results, `[[`, numeric(1), "q_long"))
  q_surv <- sum(vapply(subject_results, `[[`, numeric(1), "q_surv"))
  n_total_obs <- sum(vapply(subject_results, `[[`, numeric(1), "n_obs"))

  # Adjust longitudinal component
  q_long <- q_long - n_total_obs * half_log_2pi_sigma2

  # Random effects component (use pre-extracted values)
  inv_2tau2 <- 1 / (2 * random_effect_sd^2)
  log_2pi_tau2 <- log(2 * pi * random_effect_sd^2)
  q_re <- -sum((b_hat_all^2 + v_hat_all) * inv_2tau2) -
    length(posteriors) * log_2pi_tau2 * 0.5

  # Total Q function
  q_long + q_surv + q_re
}
