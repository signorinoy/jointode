# Utility Functions for JointODE Package
# Internal utility functions used throughout the JointODE package.

# ==============================================================================
# SECTION 1: DATA VALIDATION AND PROCESSING
# ==============================================================================

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
      if (!is.numeric(spline_baseline$boundary_knots) ||
            length(spline_baseline$boundary_knots) != 2) {
        stop(paste(
          "spline_baseline$boundary_knots must be NULL or",
          "a numeric vector of length 2"
        ))
      }
      if (spline_baseline$boundary_knots[1] >=
            spline_baseline$boundary_knots[2]) {
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
    if (!is.numeric(spline_index$degree) || length(spline_index$degree) != 1 ||
          spline_index$degree < 1 || spline_index$degree > 5) {
      stop("spline_index$degree must be a single integer between 1 and 5")
    }
  }

  if ("n_knots" %in% names(spline_index)) {
    if (!is.numeric(spline_index$n_knots) ||
          length(spline_index$n_knots) != 1 ||
          spline_index$n_knots < 1 || spline_index$n_knots > 20) {
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
      if (!is.numeric(spline_index$boundary_knots) ||
            length(spline_index$boundary_knots) != 2) {
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

# ==============================================================================
# SECTION 2: NUMERICAL COMPUTATION
# ==============================================================================

# Configure B-Spline Parameters
# Creates B-spline configuration with optimal knot placement
.get_spline_config <- function(
    x,
    degree = 3,
    n_knots = 5,
    knot_placement = "quantile",
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

# Internal function: Compute B-spline basis
.compute_spline_basis <- function(x, config) {
  splines::bs(
    x = x,
    knots = config$knots,
    degree = config$degree,
    Boundary.knots = config$boundary_knots,
    intercept = TRUE
  )
}

# Internal function: Compute B-spline basis derivative
.compute_spline_basis_deriv <- function(x, config) {
  splines::splineDesign(
    knots = c(
      rep(config$boundary_knots[1], config$degree + 1),
      config$knots,
      rep(config$boundary_knots[2], config$degree + 1)
    ),
    x = x,
    ord = config$degree + 1,
    derivs = 1
  )
}

# ==============================================================================
# SECTION 3: ODE SOLVING HELPERS
# ==============================================================================

# Internal function: Get longitudinal covariates at time point
.get_longitudinal_covariates <- function(data, time = NULL, row_index = NULL) {
  if (is.null(data$longitudinal$covariates) ||
        nrow(data$longitudinal$covariates) == 0) {
    return(numeric(0))
  }

  # If time is provided, find nearest time point
  if (!is.null(time) && !is.null(data$longitudinal$times)) {
    row_index <- which.min(abs(data$longitudinal$times - time))
  }

  # Default to last observation if no index
  if (is.null(row_index)) {
    row_index <- nrow(data$longitudinal$covariates)
  }

  as.numeric(data$longitudinal$covariates[row_index, ])
}

# Internal function: Compute acceleration function
.compute_acceleration <- function(biomarker, velocity, time_point,
                                  data, coef, config) {
  # Get longitudinal covariates
  longitudinal_covariates <- .get_longitudinal_covariates(data, time_point)

  # Build covariate vector Z
  z <- c(biomarker, velocity, longitudinal_covariates, time_point)
  n_beta <- length(coef$index_beta)
  if (length(z) < n_beta) {
    z <- c(z, rep(0, n_beta - length(z)))
  }
  z <- z[1:n_beta]

  # Compute single index
  index_value <- sum(coef$index_beta * z)

  # Evaluate acceleration function
  basis_g <- .compute_spline_basis(index_value, config$index)
  as.numeric(basis_g %*% coef$index_g)
}

# Internal function: Compute log-hazard function
.compute_log_hazard <- function(time_point, biomarker, velocity,
                                acceleration, data, coef, config) {
  # Baseline hazard
  basis_lambda <- .compute_spline_basis(time_point, config$baseline)
  log_baseline <- as.numeric(basis_lambda %*% coef$baseline)

  # Biomarker effects
  log_biomarker <- sum(
    coef$hazard[1:3] * c(biomarker, velocity, acceleration)
  )

  # Covariate effects
  log_covariate <- 0
  if (length(coef$hazard) > 3 && !is.null(data$covariates)) {
    w <- as.numeric(data$covariates)
    if (length(w) > 0) {
      log_covariate <- sum(coef$hazard[4:(3 + length(w))] * w)
    }
  }

  log_baseline + log_biomarker + log_covariate
}

# Internal function: Prepare initial conditions for ODE
.prepare_initial_conditions <- function(initial, parameters, sensitivity_type) {
  n_basis_lambda <- parameters$config$baseline$df

  switch(sensitivity_type,
    basic = initial,

    gradient = c(
      initial,
      rep(0, n_basis_lambda),
      rep(0, 3)
    ),

    beta = {
      n_beta <- length(parameters$coef$index_beta)
      c(
        initial,
        rep(0, n_basis_lambda),
        rep(0, 3),
        rep(0, 2 * n_beta),
        rep(0, n_beta)
      )
    },

    theta = {
      n_theta <- length(parameters$coef$index_g)
      c(
        initial,
        rep(0, n_basis_lambda),
        rep(0, 3),
        rep(0, 2 * n_theta),
        rep(0, n_theta)
      )
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
    solution[obs_indices, 3]  # Column 3 is m(t)
  } else {
    NULL
  }

  # Build base result
  result <- list(
    cum_hazard = final_state[1],
    log_hazard = log_hazard_final,
    biomarker = biomarker_values
  )

  # Add sensitivity-specific outputs
  if (sensitivity_type == "gradient") {
    n_basis <- parameters$config$baseline$df
    result$I_B <- final_state[4:(4 + n_basis - 1)]
    result$I_m <- final_state[(4 + n_basis):(4 + n_basis + 2)]

  } else if (sensitivity_type == "beta") {
    n_basis <- parameters$config$baseline$df
    n_beta <- length(parameters$coef$index_beta)
    idx_start <- 4 + n_basis + 3

    result$I_B <- final_state[4:(4 + n_basis - 1)]
    result$I_m <- final_state[(4 + n_basis):(4 + n_basis + 2)]
    result$S_beta <- matrix(
      final_state[(idx_start + 1):(idx_start + 2 * n_beta)],
      nrow = 2, byrow = TRUE
    )
    result$dLambda_dbeta <- final_state[
      (idx_start + 2 * n_beta + 1):(idx_start + 3 * n_beta)
    ]

  } else if (sensitivity_type == "theta") {
    n_basis <- parameters$config$baseline$df
    n_theta <- length(parameters$coef$index_g)
    idx_start <- 4 + n_basis + 3

    result$I_B <- final_state[4:(4 + n_basis - 1)]
    result$I_m <- final_state[(4 + n_basis):(4 + n_basis + 2)]
    result$S_theta <- matrix(
      final_state[(idx_start + 1):(idx_start + 2 * n_theta)],
      nrow = 2, byrow = TRUE
    )
    result$dLambda_dtheta <- final_state[
      (idx_start + 2 * n_theta + 1):(idx_start + 3 * n_theta)
    ]
  }

  result
}

# Internal function: Solve joint ODE system
# Solves the coupled ODE system for cumulative hazard and biomarker trajectory
# Supports multiple configurations: basic, gradient, beta, theta
.solve_joint_ode <- function(data, parameters, initial = c(0, 0, 0),
                             sensitivity_type = "basic") {

  # Validate sensitivity_type
  valid_types <- c("basic", "gradient", "beta", "theta")
  if (!sensitivity_type %in% valid_types) {
    stop("Invalid sensitivity_type. Must be one of: ",
         paste(valid_types, collapse = ", "))
  }

  # Define ODE derivatives function based on sensitivity type
  ode_derivatives <- if (sensitivity_type == "gradient" ||
                           sensitivity_type == "basic") {
    # ODE-1: Basic forward solve (with or without integrals)
    function(time, state, parameters) {
      biomarker_value <- state[2]
      biomarker_velocity <- state[3]

      # Get covariates at current time (nearest neighbor interpolation)
      # Note: Covariates retrieval implemented but not currently used
      # This is a placeholder for future longitudinal covariate support

      # Compute acceleration
      biomarker_acceleration <- .compute_acceleration(
        biomarker_value, biomarker_velocity, time,
        parameters$data, parameters$coef, parameters$config
      )

      # Compute hazard (without random effect b)
      log_hazard <- .compute_log_hazard(
        time, biomarker_value, biomarker_velocity, biomarker_acceleration,
        parameters$data, parameters$coef, parameters$config
      )
      hazard <- exp(log_hazard)

      # Compute B-spline basis for integrals
      basis_lambda <- .compute_spline_basis(time, parameters$config$baseline)

      # Compute m vector for integral
      m_vec <- c(biomarker_value, biomarker_velocity, biomarker_acceleration)

      # Basic derivatives: [dΛ/dt, dm/dt, dṁ/dt]
      basic_derivs <- c(hazard, biomarker_velocity, biomarker_acceleration)

      # Add integral derivatives only if sensitivity_type == "gradient"
      # (not "basic")
      if (parameters$sensitivity_type == "gradient" && length(state) > 3) {
        # I_B derivatives (positions 4 to 4+n_basis_lambda-1)
        i_b_derivs <- as.numeric(basis_lambda) * hazard

        # I_m derivatives (next 3 positions)
        i_m_derivs <- m_vec * hazard

        list(c(basic_derivs, i_b_derivs, i_m_derivs))
      } else {
        list(basic_derivs)
      }
    }

  } else if (sensitivity_type == "beta") {
    # ODE-2: With beta sensitivity
    function(time, state, parameters) {
      n_basis_lambda <- parameters$config$baseline$df
      n_beta <- length(parameters$coef$index_beta)

      # Extract base states
      biomarker_value <- state[2]
      biomarker_velocity <- state[3]

      # Extract sensitivity states
      idx_s_beta <- (4 + n_basis_lambda + 3) + 1  # Start of S_beta
      s_m_beta <- state[idx_s_beta:(idx_s_beta + n_beta - 1)]
      s_v_beta <- state[(idx_s_beta + n_beta):(idx_s_beta + 2 * n_beta - 1)]

      # Get longitudinal covariates
      has_longitudinal <- length(parameters$data$longitudinal$times) > 0
      longitudinal_covariates <- if (has_longitudinal) {
        time_index <- which.min(abs(parameters$data$longitudinal$times - time))
        .get_longitudinal_covariates(parameters$data, row_index = time_index)
      } else {
        numeric(0)
      }

      # Compute Z vector and u = beta^T Z
      z <- c(biomarker_value, biomarker_velocity, longitudinal_covariates, time)
      n_beta_actual <- length(parameters$coef$index_beta)
      if (length(z) < n_beta_actual) {
        z <- c(z, rep(0, n_beta_actual - length(z)))
      }
      z <- z[1:n_beta_actual]
      u <- sum(parameters$coef$index_beta * z)

      # Compute basis functions
      basis_g <- .compute_spline_basis(u, parameters$config$index)
      basis_g_deriv <- .compute_spline_basis_deriv(u, parameters$config$index)

      # Acceleration and its beta sensitivity
      biomarker_acceleration <- as.numeric(basis_g %*% parameters$coef$index_g)
      s_a_beta <- as.numeric(parameters$coef$index_g %*% basis_g_deriv) * z

      # Compute hazard
      log_hazard <- .compute_log_hazard(
        time, biomarker_value, biomarker_velocity, biomarker_acceleration,
        parameters$data, parameters$coef, parameters$config
      )
      hazard <- exp(log_hazard)

      # Hazard sensitivity w.r.t. beta
      dhazard_dbeta <- hazard * sum(parameters$coef$hazard[1:3] *
                                      c(s_m_beta, s_v_beta, s_a_beta))

      # Standard ODE-1 components
      basis_lambda <- .compute_spline_basis(time, parameters$config$baseline)
      m_vec <- c(biomarker_value, biomarker_velocity, biomarker_acceleration)

      # All derivatives
      basic_derivs <- c(hazard, biomarker_velocity, biomarker_acceleration)
      i_b_derivs <- as.numeric(basis_lambda) * hazard
      i_m_derivs <- m_vec * hazard
      s_beta_derivs <- c(s_v_beta, s_a_beta)
      lambda_beta_derivs <- dhazard_dbeta

      list(c(basic_derivs, i_b_derivs, i_m_derivs, s_beta_derivs,
             lambda_beta_derivs))
    }
  } else if (sensitivity_type == "theta") {
    # ODE-3: With theta sensitivity
    function(time, state, parameters) {
      n_basis_lambda <- parameters$config$baseline$df
      n_theta <- length(parameters$coef$index_g)

      # Extract base states
      biomarker_value <- state[2]
      biomarker_velocity <- state[3]

      # Extract sensitivity states
      idx_s_theta <- (4 + n_basis_lambda + 3) + 1  # Start of S_theta
      s_m_theta <- state[idx_s_theta:(idx_s_theta + n_theta - 1)]
      s_v_theta <- state[(idx_s_theta + n_theta):
                           (idx_s_theta + 2 * n_theta - 1)]

      # Get longitudinal covariates
      has_longitudinal <- length(parameters$data$longitudinal$times) > 0
      longitudinal_covariates <- if (has_longitudinal) {
        time_index <- which.min(abs(parameters$data$longitudinal$times - time))
        .get_longitudinal_covariates(parameters$data, row_index = time_index)
      } else {
        numeric(0)
      }

      # Compute u and basis
      z <- c(biomarker_value, biomarker_velocity, longitudinal_covariates, time)
      n_beta <- length(parameters$coef$index_beta)
      if (length(z) < n_beta) z <- c(z, rep(0, n_beta - length(z)))
      u <- sum(parameters$coef$index_beta * z[1:n_beta])
      basis_g <- .compute_spline_basis(u, parameters$config$index)

      # Acceleration
      biomarker_acceleration <- as.numeric(basis_g %*% parameters$coef$index_g)

      # Sensitivity of acceleration w.r.t. theta = B_g(u)
      s_a_theta <- as.numeric(basis_g)

      # Compute hazard
      log_hazard <- .compute_log_hazard(
        time, biomarker_value, biomarker_velocity, biomarker_acceleration,
        parameters$data, parameters$coef, parameters$config
      )
      hazard <- exp(log_hazard)

      # Hazard sensitivity w.r.t. theta
      dhazard_dtheta <- hazard * sum(parameters$coef$hazard[1:3] *
                                       c(s_m_theta, s_v_theta, s_a_theta))

      # Standard ODE-1 components
      basis_lambda <- .compute_spline_basis(time, parameters$config$baseline)
      m_vec <- c(biomarker_value, biomarker_velocity, biomarker_acceleration)

      # All derivatives
      basic_derivs <- c(hazard, biomarker_velocity, biomarker_acceleration)
      i_b_derivs <- as.numeric(basis_lambda) * hazard
      i_m_derivs <- m_vec * hazard
      s_theta_derivs <- c(s_v_theta, s_a_theta)
      lambda_theta_derivs <- dhazard_dtheta

      list(c(basic_derivs, i_b_derivs, i_m_derivs, s_theta_derivs,
             lambda_theta_derivs))
    }
  } else {
    stop("Unknown sensitivity type: ", sensitivity_type)
  }

  # Prepare initial conditions based on sensitivity type
  initial_extended <- .prepare_initial_conditions(initial, parameters,
                                                  sensitivity_type)

  # Solve ODE System
  ode_parameters <- list(
    data = data,
    coef = parameters$coef,
    config = parameters$config,
    sensitivity_type = sensitivity_type  # Pass sensitivity type to ODE function
  )
  event_time <- data$time
  time_points <- sort(unique(c(0, data$longitudinal$times, event_time)))

  solution <- deSolve::ode(
    y = initial_extended,
    times = time_points,
    func = ode_derivatives,
    parms = ode_parameters,
    method = "lsoda",
    rtol = 1e-8,
    atol = 1e-10
  )

  # Extract and return results
  .extract_ode_results(solution, data, parameters, sensitivity_type)
}

# ==============================================================================
# SECTION 4: STATISTICAL DISTRIBUTIONS
# ==============================================================================

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

  # Define log-posterior function (up to normalizing constant)
  logpost <- function(b) {
    ll_long <- -(n_obs * b^2 - 2 * b * s_i) / (2 * measurement_error_sd^2)
    ll_surv <- status * b - exp(b) * cum_hazard_0
    ll_prior <- -b^2 / (2 * random_effect_sd^2)
    ll_long + ll_surv + ll_prior
  }
  logpost_grad <- function(b) {
    grad_long <- (s_i - n_obs * b) / measurement_error_sd^2
    grad_surv <- status - exp(b) * cum_hazard_0
    grad_prior <- -b / random_effect_sd^2

    grad_long + grad_surv + grad_prior
  }
  logpost_hess <- function(b) {
    hess_long <- -n_obs / measurement_error_sd^2
    hess_surv <- -exp(b) * cum_hazard_0
    hess_prior <- -1 / random_effect_sd^2

    matrix(hess_long + hess_surv + hess_prior, 1, 1)
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

  # Initialize Q components
  q_long <- 0
  q_surv <- 0
  n_total_obs <- 0

  # Loop over subjects
  for (i in seq_along(data_list)) {
    subject_data <- data_list[[i]]
    subject_id <- names(data_list)[i]
    posterior <- posteriors[[subject_id]]

    # Solve ODE for current parameters
    ode_solution <- .solve_joint_ode(subject_data, parameters)

    # Extract posterior moments
    b_hat <- posterior$b_hat
    v_hat <- posterior$v_hat
    exp_b <- posterior$exp_b

    # Longitudinal component
    if (subject_data$longitudinal$n_obs > 0) {
      residuals <- subject_data$longitudinal$measurements -
        ode_solution$biomarker - b_hat
      q_long <- q_long -
        sum(residuals^2 + v_hat) / (2 * measurement_error_sd^2) -
        subject_data$longitudinal$n_obs *
          log(2 * pi * measurement_error_sd^2) / 2
      n_total_obs <- n_total_obs + subject_data$longitudinal$n_obs
    }
    q_long <- q_long - n_total_obs * log(2 * pi * measurement_error_sd^2) / 2

    # Survival component
    if (subject_data$status == 1) {
      # Event occurred
      q_surv <- q_surv + ode_solution$log_hazard + b_hat
    }
    q_surv <- q_surv - exp_b * ode_solution$cum_hazard
  }

  # Random effects component (only depends on random_effect_sd)
  q_re <- 0.0
  for (posterior in posteriors) {
    b_hat <- posterior$b_hat
    v_hat <- posterior$v_hat
    q_re <- q_re - (b_hat^2 + v_hat) / (2 * random_effect_sd^2) -
      log(2 * pi * random_effect_sd^2) / 2
  }

  # Total Q function
  q_long + q_surv + q_re
}
