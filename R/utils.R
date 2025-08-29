# Utility Functions for JointODE Package

# ===== SECTION 1: DATA VALIDATION AND PROCESSING =====
.validate <- function(
    longitudinal_formula, longitudinal_data, survival_formula, survival_data,
    id, time, spline_baseline = list(), spline_index = list()) {
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

  stopifnot(
    "Longitudinal data has no rows" = nrow(longitudinal_data) > 0,
    "Survival data has no rows" = nrow(survival_data) > 0
  )

  if (!id %in% names(longitudinal_data)) {
    stop(sprintf("ID variable '%s' not found in longitudinal data", id))
  }
  if (!time %in% names(longitudinal_data)) {
    stop(sprintf("Time variable '%s' not found in longitudinal data", time))
  }

  if (!id %in% names(survival_data)) {
    stop(sprintf("ID variable '%s' not found in survival data", id))
  }

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

  long_vars <- all.vars(longitudinal_formula)
  missing_long_vars <- setdiff(long_vars, names(longitudinal_data))
  if (length(missing_long_vars) > 0) {
    stop(sprintf(
      "Variables in longitudinal formula not found in data: %s",
      paste(missing_long_vars, collapse = ", ")
    ))
  }

  surv_vars <- all.vars(survival_formula)
  missing_surv_vars <- setdiff(surv_vars, names(survival_data))
  if (length(missing_surv_vars) > 0) {
    stop(sprintf(
      "Variables in survival formula not found in data: %s",
      paste(missing_surv_vars, collapse = ", ")
    ))
  }


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

  if (any(duplicated(survival_data[[id]]))) {
    stop(paste(
      "Duplicate IDs found in survival data -",
      "each subject should have one record"
    ))
  }

  long_ids <- unique(longitudinal_data[[id]])
  surv_ids <- unique(survival_data[[id]])

  orphaned_ids <- setdiff(long_ids, surv_ids)
  if (length(orphaned_ids) > 0) {
    stop(sprintf(
      "Subjects in longitudinal data not found in survival data: %s",
      paste(head(orphaned_ids, 5), collapse = ", ")
    ))
  }

  missing_long <- setdiff(surv_ids, long_ids)
  if (length(missing_long) > 0) {
    warning(sprintf(
      "Subjects in survival data without longitudinal data: %s",
      paste(head(missing_long, 5), collapse = ", ")
    ))
  }


  if (any(longitudinal_data[[time]] < 0, na.rm = TRUE)) {
    stop("Negative time values found in longitudinal data")
  }

  if (any(survival_data[[time_var]] <= 0, na.rm = TRUE)) {
    stop("Invalid observation times in survival data (must be positive)")
  }

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

  observations_per_subject <- table(longitudinal_data[[id]])
  if (all(observations_per_subject == 1)) {
    warning(paste(
      "Each subject has only one longitudinal observation -",
      "joint modeling may not be appropriate"
    ))
  }

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

.process <- function(
    longitudinal_formula, longitudinal_data, survival_formula, survival_data,
    id, time) {
  unique_ids <- unique(survival_data[[id]])
  n_subjects <- length(unique_ids)
  data_process <- vector("list", n_subjects)
  names(data_process) <- as.character(unique_ids)

  surv_frame <- model.frame(survival_formula, data = survival_data)
  surv_response_matrix <- model.response(surv_frame)
  has_surv_covs <- length(all.vars(survival_formula[[3]])) > 0 &&
    survival_formula[[3]] != 1
  surv_design <- if (has_surv_covs) {
    model.matrix(survival_formula, surv_frame)[, -1, drop = FALSE]
  } else {
    NULL
  }


  survival_index_map <- match(unique_ids, survival_data[[id]])
  long_id_groups <- split(
    seq_len(nrow(longitudinal_data)), longitudinal_data[[id]]
  )

  for (i in seq_along(unique_ids)) {
    sid <- as.character(unique_ids[i])

    long_rows <- long_id_groups[[sid]]
    if (!is.null(long_rows) && length(long_rows) > 0) {
      long_subset <- longitudinal_data[long_rows, , drop = FALSE]
      ord <- order(long_subset[[time]])
      if (!identical(ord, seq_along(ord))) {
        long_subset <- long_subset[ord, , drop = FALSE]
      }
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

    survival_row <- survival_index_map[i]
    event_time <- surv_response_matrix[survival_row, 1]
    event_status <- surv_response_matrix[survival_row, 2]
    covariates <- if (!is.null(surv_design)) {
      surv_design[survival_row, , drop = FALSE]
    } else {
      data.frame()
    }

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

  statuses <- vapply(data_process, `[[`, numeric(1), "status")
  attr(data_process, "n_subjects") <- n_subjects
  attr(data_process, "n_observations") <- nrow(longitudinal_data)
  attr(data_process, "event_rate") <- mean(statuses == 1, na.rm = TRUE)

  data_process
}

.get_spline_config <- function(
    x, degree = 3, n_knots = 5, knot_placement = "quantile",
    boundary_knots = NULL) {
  if (is.null(boundary_knots)) {
    boundary_knots <- range(x, na.rm = TRUE)
  } else {
    boundary_knots <- boundary_knots
  }

  if (knot_placement == "quantile") {
    probs <- seq(0, 1, length.out = n_knots + 2)[-c(1, n_knots + 2)]
    knots <- quantile(x, probs = probs, na.rm = TRUE, names = FALSE)
  } else if (knot_placement == "equal") {
    knots <- seq(boundary_knots[1], boundary_knots[2],
      length.out = n_knots + 2
    )[-c(1, n_knots + 2)]
  } else {
    stop("knot_placement must be 'quantile' or 'equal'")
  }

  list(
    degree = degree, knots = knots, boundary_knots = boundary_knots,
    df = length(knots) + degree + 1
  )
}

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

.compute_spline_basis_deriv <- function(x, config) {
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


.get_longitudinal_covariates <- function(data, time = NULL, row_index = NULL) {
  long_cov <- data$longitudinal$covariates
  if (is.null(long_cov) || nrow(long_cov) == 0) {
    return(numeric(0))
  }

  if (is.null(row_index)) {
    if (!is.null(time) && !is.null(data$longitudinal$times)) {
      row_index <- which.min(abs(data$longitudinal$times - time))
    } else {
      row_index <- nrow(long_cov)
    }
  }

  as.numeric(long_cov[row_index, , drop = TRUE])
}

.compute_acceleration <- function(
    biomarker, velocity, time, data, parameters) {
  long_cov <- .get_longitudinal_covariates(data, time)
  z <- as.vector(c(biomarker, velocity, long_cov, time))
  index_value <- sum(parameters$coefficients$index_beta * z)

  basis_g <- .compute_spline_basis(index_value, parameters$configurations$index)
  sum(basis_g * parameters$coefficients$index_g)
}

.compute_acceleration_deriv <- function(
    biomarker, velocity, time, data, parameters, type = "beta",
    dbiomarker = NULL, dvelocity = NULL) {
  long_cov <- .get_longitudinal_covariates(data, time)

  # Construct Z vector: [m(t), m_dot(t), X(t), t]
  z_vec <- c(biomarker, velocity, long_cov, time)

  beta <- parameters$coefficients$index_beta
  theta <- parameters$coefficients$index_g

  # Compute index value u = beta' * Z
  index_value <- sum(beta * z_vec)

  # Compute B-spline basis and its derivative
  basis_g <- .compute_spline_basis(
    index_value, parameters$configurations$index
  )
  basis_g_deriv <- .compute_spline_basis_deriv(
    index_value, parameters$configurations$index
  )

  if (type == "beta") {
    # Compute d(acceleration)/dbeta
    # d(ddot{m})/dbeta = theta' * B'_g(u) * du/dbeta
    # where du/dbeta = Z + beta' * dZ/dbeta

    n_beta <- length(beta)

    if (!is.null(dbiomarker) && !is.null(dvelocity)) {
      du_dbeta <- numeric(n_beta)
      for (j in seq_len(n_beta)) {
        # Direct effect: partial u / partial beta_j
        direct_effect <- if (j <= length(z_vec)) z_vec[j] else 0
        # Indirect effect through Z: beta1 * dm/dbeta_j + beta2 * dm'/dbeta_j
        indirect_effect <- beta[1] * dbiomarker[j] + beta[2] * dvelocity[j]
        du_dbeta[j] <- direct_effect + indirect_effect
      }
    } else {
      # If sensitivities not provided, use just the direct effect
      du_dbeta <- numeric(n_beta)
      for (j in seq_len(n_beta)) {
        du_dbeta[j] <- if (j <= length(z_vec)) z_vec[j] else 0
      }
    }

    # Compute d(acceleration)/dbeta = theta' * B'_g(u) * du/dbeta
    dacceleration_dbeta <- sum(theta * basis_g_deriv) * du_dbeta

    return(dacceleration_dbeta)
  } else if (type == "theta") {
    # Compute d(acceleration)/dgamma for gamma parameters
    # (theta in this context)
    # acceleration = gamma' * B_g(u) where u = beta' * Z
    # For each gamma_j:
    # d(acceleration)/dgamma_j = B_g_j(u) + gamma' * B'_g(u) * du/dgamma_j
    # where du/dgamma_j = beta[1] * dm/dgamma_j + beta[2] * dṁ/dgamma_j

    # Compute the coefficient for the indirect effect
    indirect_coef <- sum(theta * basis_g_deriv)

    # Initialize result vector
    n_gamma <- length(dbiomarker)
    dacceleration_dgamma <- numeric(n_gamma)

    # For each gamma component
    for (j in 1:n_gamma) {
      # du/dgamma_j for this specific gamma_j
      du_dgamma_j <- beta[1] * dbiomarker[j] + beta[2] * dvelocity[j]

      # Total derivative for gamma_j
      # Direct effect: B_g_j(u) (only non-zero for component j)
      # Indirect effect: gamma' * B'_g(u) * du/dgamma_j
      dacceleration_dgamma[j] <- basis_g[j] + indirect_coef * du_dgamma_j
    }

    return(dacceleration_dgamma)
  } else {
    stop("Invalid type. Must be one of: 'beta', 'theta'")
  }
}

.compute_log_hazard <- function(
    time, biomarker, velocity, acceleration, data, parameters) {
  # Baseline hazard
  basis_lambda <- .compute_spline_basis(
    time, parameters$configurations$baseline
  )
  log_baseline <- sum(basis_lambda * parameters$coefficients$baseline)

  # Biomarker effects: m_i(t)' * α
  biomarker_vec <- c(biomarker, velocity, acceleration)
  log_biomarker <- sum(biomarker_vec * parameters$coefficients$hazard[1:3])

  # Covariate effects: W_i' * φ
  n_hazard <- length(parameters$coefficients$hazard)
  log_covariate <- if (n_hazard > 3 && !is.null(data$covariates)) {
    w <- data$covariates
    if (nrow(w) > 0) {
      sum(parameters$coefficients$hazard[4:n_hazard] * w)
    } else {
      0
    }
  } else {
    0
  }
  log_baseline + log_biomarker + log_covariate
}

.prepare_initial_conditions <- function(
    biomarker_initial, coefficients, sensitivity_type) {
  if (length(biomarker_initial) != 2) {
    stop("biomarker_initial must be a numeric vector of length 2")
  }

  # Basic state: [Λ(t), m(t), ṁ(t)]
  basic_state <- c(0, biomarker_initial)

  switch(sensitivity_type,
    basic = basic_state,
    forward_beta = {
      # Augmented for β sensitivities:
      # [Λ, m, ṁ, ∂m/∂β, ∂ṁ/∂β, ∂Λ/∂β]
      n_beta <- length(coefficients$coef$index_beta)
      c(
        basic_state,
        rep(0, n_beta), # ∂m/∂β
        rep(0, n_beta), # ∂ṁ/∂β
        rep(0, n_beta)  # ∂Λ/∂β
      )
    },
    forward_theta = {
      # Combined augmented system for θ = (η, α, γ) sensitivities:
      # [Λ, m, ṁ, ∂Λ/∂η, ∂Λ/∂α, ∂m/∂γ, ∂ṁ/∂γ, ∂Λ/∂γ]
      n_eta <- coefficients$config$baseline$df
      n_alpha <- 3 # [α₀, α₁, α₂] for biomarker, velocity, acceleration
      n_gamma <- length(coefficients$coef$index_g)
      c(
        basic_state,
        rep(0, n_eta),   # ∂Λ/∂η
        rep(0, n_alpha), # ∂Λ/∂α
        rep(0, n_gamma), # ∂m/∂γ
        rep(0, n_gamma), # ∂ṁ/∂γ
        rep(0, n_gamma)  # ∂Λ/∂γ
      )
    },
    adjoint_beta = {
      # For adjoint method (future implementation)
      # Initial conditions will be set at terminal time
      basic_state
    },
    adjoint_theta = {
      # For adjoint method (future implementation)
      # Initial conditions will be set at terminal time
      basic_state
    },
    {
      stop("Invalid sensitivity_type: ", sensitivity_type)
    }
  )
}

.extract_ode_results <- function(solution, data, parameters, sensitivity_type) {
  # Get final state
  final_state <- solution[nrow(solution), -1]
  event_time <- data$time

  # Compute final values
  biomarker_final <- final_state[2]
  velocity_final <- final_state[3]

  acceleration_final <- .compute_acceleration(
    biomarker_final, velocity_final, event_time, data, parameters
  )

  log_hazard_final <- .compute_log_hazard(
    event_time, biomarker_final, velocity_final, acceleration_final,
    data, parameters
  )

  # Extract biomarker trajectory at observation times
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
        solution[idx, 3], solution[idx, 4], solution[idx, 1], data, parameters
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
  if (sensitivity_type == "forward_theta") {
    # Combined θ = (η, α, γ) sensitivities
    n_eta <- parameters$configurations$baseline$df
    n_alpha <- 3
    n_gamma <- length(parameters$coefficients$index_g)

    # Extract sensitivities from final state (matching ODE output order)
    # State order: [Λ, m, ṁ, ∂Λ/∂η, ∂Λ/∂α, ∂m/∂γ, ∂ṁ/∂γ, ∂Λ/∂γ]
    idx <- 4
    dcumhazard_deta <- final_state[idx:(idx + n_eta - 1)]
    idx <- idx + n_eta
    dcumhazard_dalpha <- final_state[idx:(idx + n_alpha - 1)]
    idx <- idx + n_alpha
    dbiomarker_dgamma <- final_state[idx:(idx + n_gamma - 1)]
    idx <- idx + n_gamma
    dvelocity_dgamma <- final_state[idx:(idx + n_gamma - 1)]
    idx <- idx + n_gamma
    dcumhazard_dgamma <- final_state[idx:(idx + n_gamma - 1)]

    # Extract γ sensitivities at observation times
    dbiomarker_dgamma_values <- if (length(data$longitudinal$times) > 0) {
      obs_indices <- match(data$longitudinal$times, solution[, 1])
      # Column 1 is time, columns 2-4 are [Λ, m, ṁ]
      # Columns 5-(4+n_eta) are ∂Λ/∂η
      # Columns (5+n_eta)-(4+n_eta+n_alpha) are ∂Λ/∂α
      # Columns (5+n_eta+n_alpha)-(4+n_eta+n_alpha+n_gamma) are ∂m/∂γ
      col_start <- 5 + n_eta + n_alpha  # Adjust for time column (column 1)
      col_end <- col_start + n_gamma - 1
      solution[obs_indices, col_start:col_end, drop = FALSE]
    } else {
      NULL
    }

    # Compute acceleration sensitivity at event time
    dacceleration_dgamma <- .compute_acceleration_deriv(
      biomarker_final, velocity_final, event_time, data, parameters,
      type = "theta",
      dbiomarker = dbiomarker_dgamma,
      dvelocity = dvelocity_dgamma
    )

    result <- list(
      log_hazard = log_hazard_final,
      cum_hazard = final_state[1],
      biomarker_at_event = final_state[2],
      velocity_at_event = final_state[3],
      acceleration_at_event = acceleration_final,
      biomarker = biomarker_values,
      # η sensitivities
      dcumhazard_deta_at_event = dcumhazard_deta,
      # α sensitivities
      dcumhazard_dalpha_at_event = dcumhazard_dalpha,
      # γ sensitivities
      dbiomarker_dgamma = dbiomarker_dgamma_values,
      dbiomarker_dgamma_at_event = dbiomarker_dgamma,
      dvelocity_dgamma_at_event = dvelocity_dgamma,
      dacceleration_dgamma_at_event = dacceleration_dgamma,
      dcumhazard_dgamma_at_event = dcumhazard_dgamma
    )
  } else if (sensitivity_type == "forward_beta") {
    # β sensitivities
    n_beta <- length(parameters$coefficients$index_beta)

    # Extract β sensitivities at observation times
    dbiomarker_dbeta_values <- if (length(data$longitudinal$times) > 0) {
      obs_indices <- match(data$longitudinal$times, solution[, 1])
      solution[obs_indices, 5:(4 + n_beta), drop = FALSE]
    } else {
      NULL
    }

    # Extract final sensitivities
    dbiomarker_dbeta_final <- final_state[4:(3 + n_beta)]
    dvelocity_dbeta_final <- final_state[(4 + n_beta):(3 + 2 * n_beta)]
    dcumhazard_dbeta_final <- final_state[(4 + 2 * n_beta):(3 + 3 * n_beta)]

    # Compute acceleration sensitivity at event time
    dacceleration_dbeta_final <- .compute_acceleration_deriv(
      biomarker_final, velocity_final, event_time, data, parameters,
      type = "beta",
      dbiomarker = dbiomarker_dbeta_final,
      dvelocity = dvelocity_dbeta_final
    )

    result <- list(
      log_hazard = log_hazard_final,
      cum_hazard = final_state[1],
      biomarker = biomarker_values,
      dbiomarker_dbeta = dbiomarker_dbeta_values,
      dbiomarker_dbeta_at_event = dbiomarker_dbeta_final,
      dvelocity_dbeta_at_event = dvelocity_dbeta_final,
      dacceleration_dbeta_at_event = dacceleration_dbeta_final,
      dcumhazard_dbeta_at_event = dcumhazard_dbeta_final
    )
  }

  result
}

.solve_joint_ode <- function(data, parameters, sensitivity_type = "basic") {
  # Validate sensitivity_type
  valid_types <- c(
    "basic",
    "forward_beta", "forward_theta",
    "adjoint_beta", "adjoint_theta"
  )
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
      biomarker, velocity, time, parameters$data, parameters$parameters
    )

    # Compute hazard (without random effect b)
    log_hazard <- .compute_log_hazard(
      time, biomarker, velocity, acceleration,
      parameters$data, parameters$parameters
    )
    # Cap log_hazard to prevent numerical overflow
    log_hazard_capped <- pmin(log_hazard, 50) # exp(50) ≈ 5e21
    hazard <- exp(log_hazard_capped)

    # Basic derivatives: [dΛ/dt, dm/dt, dṁ/dt]
    basic_derivs <- c(hazard, velocity, acceleration)

    if (parameters$sensitivity_type == "basic") {
      list(basic_derivs)
    } else if (parameters$sensitivity_type == "forward_theta") {
      # Combined forward sensitivity for θ = (η, α, γ)
      n_eta <- parameters$parameters$configurations$baseline$df
      n_alpha <- 3
      n_gamma <- length(parameters$parameters$coefficients$index_g)

      # Compute ∂Λ/∂η and ∂Λ/∂α (direct effects on hazard)
      basis_lambda <- .compute_spline_basis(
        time, parameters$parameters$configurations$baseline
      )
      m_vec <- c(biomarker, velocity, acceleration)

      # Extract γ sensitivities from state
      idx <- 4 + n_eta + n_alpha
      dbiomarker_dgamma <- state[idx:(idx + n_gamma - 1)]
      idx <- idx + n_gamma
      dvelocity_dgamma <- state[idx:(idx + n_gamma - 1)]

      # Compute acceleration sensitivity w.r.t. γ
      dacceleration_dgamma <- .compute_acceleration_deriv(
        biomarker, velocity, time, parameters$data, parameters$parameters,
        type = "theta",
        dbiomarker = dbiomarker_dgamma,
        dvelocity = dvelocity_dgamma
      )

      # Compute ∂λ/∂γ (hazard sensitivity w.r.t. gamma) using chain rule
      alpha <- parameters$parameters$coefficients$hazard[1:3]
      dm_vec_dgamma <- rbind(
        dbiomarker_dgamma,
        dvelocity_dgamma,
        dacceleration_dgamma
      )
      dhazard_dgamma <- as.vector(
        crossprod(alpha, dm_vec_dgamma)
      ) * hazard

      # Return augmented derivatives
      # State order: [Λ, m, ṁ, ∂Λ/∂η, ∂Λ/∂α, ∂m/∂γ, ∂ṁ/∂γ, ∂Λ/∂γ]
      # Derivative order: [dΛ/dt, dm/dt, dṁ/dt, d(∂Λ/∂η)/dt, d(∂Λ/∂α)/dt,
      #                    d(∂m/∂γ)/dt, d(∂ṁ/∂γ)/dt, d(∂Λ/∂γ)/dt]
      list(c(
        basic_derivs, # [dΛ/dt, dm/dt, dṁ/dt]
        as.numeric(basis_lambda) * hazard, # d(∂Λ/∂η)/dt = ∂λ/∂η
        m_vec * hazard, # d(∂Λ/∂α)/dt = ∂λ/∂α
        dvelocity_dgamma, # d(∂m/∂γ)/dt = ∂ṁ/∂γ (since dm/dt = ṁ)
        dacceleration_dgamma, # d(∂ṁ/∂γ)/dt = ∂m̈/∂γ (since dṁ/dt = m̈)
        dhazard_dgamma # d(∂Λ/∂γ)/dt = ∂λ/∂γ
      ))
    } else if (parameters$sensitivity_type == "forward_beta") {
      # Forward sensitivity for single-index coefficients (β)
      n_beta <- length(parameters$parameters$coefficients$index_beta)

      # Extract sensitivity states
      dbiomarker_dbeta <- state[4:(3 + n_beta)]
      dvelocity_dbeta <- state[(4 + n_beta):(3 + 2 * n_beta)]

      # Compute acceleration sensitivity
      dacceleration_dbeta <- .compute_acceleration_deriv(
        biomarker, velocity, time, parameters$data, parameters$parameters,
        type = "beta",
        dbiomarker = dbiomarker_dbeta,
        dvelocity = dvelocity_dbeta
      )

      # Compute ∂λ/∂β (hazard sensitivity) using chain rule
      alpha <- parameters$parameters$coefficients$hazard[1:3]
      dm_vec_dbeta <- rbind(
        dbiomarker_dbeta,
        dvelocity_dbeta,
        dacceleration_dbeta
      )
      dhazard_dbeta <- as.vector(
        crossprod(alpha, dm_vec_dbeta)
      ) * hazard

      # Return augmented derivatives: [basic; ∂ṁ/∂β; ∂m̈/∂β; ∂Λ/∂β]
      list(c(
        basic_derivs, # [dΛ/dt, dm/dt, dṁ/dt]
        dvelocity_dbeta, # ∂ṁ/∂β
        dacceleration_dbeta, # ∂m̈/∂β
        dhazard_dbeta # ∂λ/∂β
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
    parameters = parameters,
    sensitivity_type = sensitivity_type
  )
  event_time <- data$time
  times <- sort(unique(c(0, data$longitudinal$times, event_time)))

  # Solve ODE with robust error handling and higher precision
  solution <- deSolve::ode(
    y = initial_extended,
    times = times,
    func = ode_derivatives,
    parms = ode_parameters,
    method = "ode45",
    atol = 1e-10,  # Increased precision (was 1e-3)
    rtol = 1e-12   # Increased precision (was 1e-4)
  )

  # Extract and return results
  .extract_ode_results(solution, data, parameters, sensitivity_type)
}

# Helper Functions for Objective and Gradient Computation
.parse_theta_parameters <- function(params, configurations, n_surv_covariates) {
  # Parse θ = (η, α, φ, γ) parameters
  n_eta <- configurations$baseline$df
  n_alpha <- 3
  n_phi <- n_surv_covariates
  n_gamma <- configurations$index$df

  idx <- 1
  eta <- params[idx:(idx + n_eta - 1)]
  idx <- idx + n_eta
  alpha <- params[idx:(idx + n_alpha - 1)]
  idx <- idx + n_alpha
  phi <- if (n_phi > 0) params[idx:(idx + n_phi - 1)] else numeric(0)
  idx <- idx + n_phi
  gamma <- params[idx:(idx + n_gamma - 1)]

  list(
    eta = eta, alpha = alpha, phi = phi, gamma = gamma,
    n_eta = n_eta, n_alpha = n_alpha, n_phi = n_phi, n_gamma = n_gamma
  )
}

# Objective Function and Gradient Computation

.compute_objective_theta <- function(
    params, data_list, posteriors, configurations, fixed_params) {
  # Compute objective function for θ = (η, α, φ, γ)
  # Returns negative expected complete-data log-likelihood (for minimization)

  n_subjects <- length(data_list)
  n_surv_covariates <- ncol(data_list[[1]]$covariates)

  # Parse parameters using helper
  theta_params <- .parse_theta_parameters(params, configurations, n_surv_covariates)

  # Build parameters structure
  parameters <- list(
    coefficients = list(
      baseline = theta_params$eta,
      hazard = c(theta_params$alpha, theta_params$phi),
      index_g = theta_params$gamma,
      index_beta = fixed_params$index_beta
    ),
    configurations = configurations
  )

  # Measurement error variance
  sigma_e <- fixed_params$measurement_error_sd
  inv_sigma_e2 <- 1 / (sigma_e^2)

  # Compute Q value
  q_value <- 0

  for (i in seq_len(n_subjects)) {
    subject_data <- data_list[[i]]

    # Solve basic ODE (no sensitivities needed)
    ode_sol <- .solve_joint_ode(
      subject_data, parameters,
      sensitivity_type = "basic"
    )

    # Survival component
    q_value <- q_value + data_list[[i]]$status * ode_sol$log_hazard -
      posteriors$exp_b[i] * ode_sol$cum_hazard

    # Longitudinal component
    if (subject_data$longitudinal$n_obs > 0) {
      residuals <- subject_data$longitudinal$measurements -
        ode_sol$biomarker - posteriors$b[i]
      q_value <- q_value - sum(residuals^2) * inv_sigma_e2 / 2
    }
  }

  # Return negative Q for minimization
  -q_value
}

.compute_objective_beta <- function(
    params, data_list, posteriors, configurations, fixed_params) {
  # Compute objective function for β (with spherical parameterization)

  # Convert spherical coordinates to unit vector
  beta <- .spherical_to_beta(params)

  # Get number of subjects
  n_subjects <- length(data_list)

  # Build parameters structure
  parameters <- list(
    coefficients = list(
      baseline = fixed_params$baseline,
      hazard = fixed_params$hazard,
      index_g = fixed_params$index_g,
      index_beta = beta
    ),
    configurations = configurations
  )

  # Measurement error variance
  sigma_e <- fixed_params$measurement_error_sd
  inv_sigma_e2 <- 1 / (sigma_e^2)

  # Compute Q value
  q_value <- 0

  for (i in seq_len(n_subjects)) {
    subject_data <- data_list[[i]]

    # Solve basic ODE
    ode_sol <- .solve_joint_ode(
      subject_data, parameters,
      sensitivity_type = "basic"
    )

    # Longitudinal component
    if (subject_data$longitudinal$n_obs > 0) {
      residuals <- subject_data$longitudinal$measurements -
        ode_sol$biomarker - posteriors$b[i]
      q_value <- q_value - sum(residuals^2) * inv_sigma_e2 / 2
    }

    # Survival component
    q_value <- q_value + data_list[[i]]$status * ode_sol$log_hazard -
      posteriors$exp_b[i] * ode_sol$cum_hazard
  }

  # Return negative Q for minimization
  -q_value
}

.compute_grad_theta_forward <- function(
    params, data_list, posteriors, configurations, fixed_params) {
  # Gradient computation for θ = (η, α, φ, γ) using forward sensitivity

  n_subjects <- length(data_list)
  n_surv_covariates <- ncol(data_list[[1]]$covariates)

  # Parse parameters using helper
  theta_params <- .parse_theta_parameters(
    params, configurations, n_surv_covariates
  )

  # Build parameters structure
  parameters <- list(
    coefficients = list(
      baseline = theta_params$eta,
      hazard = c(theta_params$alpha, theta_params$phi),
      index_g = theta_params$gamma,
      index_beta = fixed_params$index_beta
    ),
    configurations = configurations
  )

  # Initialize gradient components
  grad_eta <- numeric(theta_params$n_eta)
  grad_alpha <- numeric(theta_params$n_alpha)
  grad_phi <- numeric(theta_params$n_phi)
  grad_gamma <- numeric(theta_params$n_gamma)

  # Measurement error variance
  sigma_e <- fixed_params$measurement_error_sd
  inv_sigma_e2 <- 1 / (sigma_e^2)

  # Solve ODE with forward sensitivities for each subject
  for (i in seq_len(n_subjects)) {
    subject_data <- data_list[[i]]

    # Solve augmented ODE system
    ode_sol <- .solve_joint_ode(
      subject_data, parameters,
      sensitivity_type = "forward_theta"
    )

    status_i <- subject_data$status
    exp_b_i <- posteriors$exp_b[i]
    b_hat_i <- posteriors$b[i]

    # Longitudinal gradient (γ only)
    if (subject_data$longitudinal$n_obs > 0) {
      residuals <- subject_data$longitudinal$measurements -
        ode_sol$biomarker - b_hat_i

      if (theta_params$n_gamma > 0 && !is.null(ode_sol$dbiomarker_dgamma)) {
        grad_gamma <- grad_gamma +
          as.vector(crossprod(residuals, ode_sol$dbiomarker_dgamma)) *
            inv_sigma_e2
      }
    }

    # η gradient (baseline hazard)
    if (status_i == 1) {
      basis_at_event <- .compute_spline_basis(
        subject_data$time, configurations$baseline
      )
      grad_eta <- grad_eta + as.vector(basis_at_event)
    }
    grad_eta <- grad_eta - exp_b_i * ode_sol$dcumhazard_deta_at_event

    # α gradient (association parameters)
    if (status_i == 1) {
      m_vec_at_event <- c(
        ode_sol$biomarker_at_event,
        ode_sol$velocity_at_event,
        ode_sol$acceleration_at_event
      )
      grad_alpha <- grad_alpha + m_vec_at_event
    }
    grad_alpha <- grad_alpha - exp_b_i * ode_sol$dcumhazard_dalpha_at_event

    # φ gradient (survival covariates)
    if (theta_params$n_phi > 0 && !is.null(subject_data$covariates)) {
      covs <- as.vector(subject_data$covariates)
      grad_phi <- grad_phi + (status_i - exp_b_i * ode_sol$cum_hazard) * covs
    }

    # γ gradient (acceleration spline) - survival contribution
    if (theta_params$n_gamma > 0) {
      if (status_i == 1) {
        grad_gamma <- grad_gamma +
          theta_params$alpha[1] * ode_sol$dbiomarker_dgamma_at_event +
          theta_params$alpha[2] * ode_sol$dvelocity_dgamma_at_event +
          theta_params$alpha[3] * ode_sol$dacceleration_dgamma_at_event
      }
      grad_gamma <- grad_gamma - exp_b_i * ode_sol$dcumhazard_dgamma_at_event
    }
  }

  # Return negative gradient for minimization
  grad_full <- -c(grad_eta, grad_alpha, grad_phi, grad_gamma)
  names(grad_full) <- NULL  # Remove names to match expected format
  grad_full
}

.compute_grad_beta_forward <- function(
    params, data_list, posteriors, configurations, fixed_params) {
  # Gradient computation for β using forward sensitivity
  # params are spherical coordinates, we convert to unit vector

  beta <- .spherical_to_beta(params)
  n_beta <- length(beta)
  n_subjects <- length(data_list)

  # Extract fixed parameters
  sigma_e <- fixed_params$measurement_error_sd
  inv_sigma_e2 <- 1 / (sigma_e^2)
  alpha <- fixed_params$hazard[1:3]

  # Build parameters structure
  parameters <- list(
    coefficients = list(
      baseline = fixed_params$baseline,
      hazard = fixed_params$hazard,
      index_g = fixed_params$index_g,
      index_beta = beta
    ),
    configurations = configurations
  )
  # Initialize gradient
  grad_beta <- numeric(n_beta)

  # Process each subject
  for (i in seq_len(n_subjects)) {
    subject_data <- data_list[[i]]

    # Solve ODE with β sensitivities
    ode_sol <- .solve_joint_ode(
      subject_data, parameters,
      sensitivity_type = "forward_beta"
    )

    # Longitudinal gradient
    if (subject_data$longitudinal$n_obs > 0) {
      residuals <- subject_data$longitudinal$measurements -
        ode_sol$biomarker - posteriors$b[i]
      grad_beta <- grad_beta +
        as.vector(crossprod(residuals, ode_sol$dbiomarker_dbeta)) * inv_sigma_e2
    }

    # Survival gradient
    if (data_list[[i]]$status == 1) {
      grad_beta <- grad_beta +
        alpha[1] * ode_sol$dbiomarker_dbeta_at_event +
        alpha[2] * ode_sol$dvelocity_dbeta_at_event +
        alpha[3] * ode_sol$dacceleration_dbeta_at_event
    }
    grad_beta <- grad_beta -
      posteriors$exp_b[i] * ode_sol$dcumhazard_dbeta_at_event
  }

  # Transform gradient from β space to spherical coordinate space
  J <- .spherical_jacobian(params)
  grad_spherical <- as.vector(crossprod(J, -grad_beta))

  # Return negative gradient for minimization (in spherical coordinates)
  grad_spherical
}


# Statistical Distributions

.compute_posteriors <- function(data_list, parameters, k = 7, init = NULL) {
  n_subjects <- length(data_list)
  if (is.null(init)) init <- rep(0, n_subjects)

  posteriors <- list(
    b = numeric(n_subjects),
    v = numeric(n_subjects),
    exp_b = numeric(n_subjects)
  )

  for (i in 1:n_subjects) {
    ode_solution <- .solve_joint_ode(data_list[[i]], parameters)
    posterior <- .compute_posterior_aghq(
      ode_solution = ode_solution,
      data = data_list[[i]],
      b_hat_init = init[i],
      measurement_error_sd = parameters$coefficients$measurement_error_sd,
      random_effect_sd = parameters$coefficients$random_effect_sd,
      k = k
    )
    posteriors$b[i] <- posterior$b
    posteriors$v[i] <- posterior$v
    posteriors$exp_b[i] <- posterior$exp_b
  }

  posteriors
}

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
    b = b_hat,
    v = v_hat,
    exp_b = exp_b
  )
}

.compute_sds <- function(data_list, parameters, posteriors) {
  n_subjects <- length(data_list)
  n_observations <- sum(sapply(data_list, function(d) d$longitudinal$n_obs))
  measurement_error_variance <- 0
  random_effect_variance <- 0
  for (i in seq_along(data_list)) {
    n_i <- data_list[[i]]$longitudinal$n_obs
    ode_sol <- .solve_joint_ode(data_list[[i]], parameters)
    if (n_i > 0) {
      residuals <- data_list[[i]]$longitudinal$measurements -
        posteriors$b[i] -  ode_sol$biomarker
      measurement_error_variance <- measurement_error_variance +
        sum(residuals^2) + n_i * posteriors$v[i]
      random_effect_variance <- random_effect_variance +
        posteriors$b[i]^2 + posteriors$v[i]
    } else {
      random_effect_variance <- random_effect_variance +
        posteriors$b[i]^2 + posteriors$v[i]
    }
  }
  measurement_error_sd <- sqrt(measurement_error_variance / n_observations)
  random_effect_sd <- sqrt(random_effect_variance / n_subjects)

  list(
    measurement_error_sd = measurement_error_sd,
    random_effect_sd = random_effect_sd
  )
}

# Convert spherical coordinates to unit vector
# theta: vector of n-1 angles for n-dimensional unit vector
.spherical_to_beta <- function(theta) {
  n <- length(theta) + 1
  beta <- numeric(n)

  if (n == 1) {
    beta[1] <- 1
    return(beta)
  }

  # Convert spherical coordinates to Cartesian
  beta[1] <- cos(theta[1])

  if (n > 2) {
    for (i in 2:(n - 1)) {
      beta[i] <- prod(sin(theta[1:(i - 1)])) * cos(theta[i])
    }
  }

  beta[n] <- prod(sin(theta))

  beta
}

# Convert unit vector to spherical coordinates
.beta_to_spherical <- function(beta) {
  n <- length(beta)
  if (n == 1) {
    return(numeric(0))
  }

  theta <- numeric(n - 1)

  # First angle
  theta[1] <- acos(pmin(pmax(beta[1], -1), 1))

  if (n > 2) {
    for (i in 2:(n - 1)) {
      sin_prod <- prod(sin(theta[1:(i - 1)]))
      if (abs(sin_prod) > 1e-10) {
        theta[i] <- acos(pmin(pmax(beta[i] / sin_prod, -1), 1))
      } else {
        theta[i] <- 0
      }
    }
  }

  # Adjust last angle based on sign of last component
  if (beta[n] < 0) {
    theta[n - 1] <- 2 * pi - theta[n - 1]
  }

  theta
}

# Jacobian matrix: d(beta)/d(theta)
.spherical_jacobian <- function(theta) {
  n <- length(theta) + 1
  J <- matrix(0, n, length(theta))

  if (n == 1) {
    return(J)
  }

  # First component
  J[1, 1] <- -sin(theta[1])

  if (n > 2) {
    for (i in 2:(n - 1)) {
      # Component i depends on theta[1] through theta[i]
      for (j in seq_len(min(i, length(theta)))) {
        if (j < i) {
          # Derivative w.r.t. theta[j] for j < i
          if (j == 1) {
            sin_prod <- if (i > 2) prod(sin(theta[2:(i - 1)])) else 1
            J[i, j] <- cos(theta[1]) * sin_prod * cos(theta[i])
          } else {
            prod_before <- if (j > 2) {
              prod(sin(theta[1:(j - 1)]))
            } else {
              sin(theta[1])
            }
            prod_after <- if (j < i - 1) {
              prod(sin(theta[(j + 1):(i - 1)]))
            } else {
              1
            }
            J[i, j] <- prod_before * cos(theta[j]) * prod_after * cos(theta[i])
          }
        } else if (j == i) {
          # Derivative w.r.t. theta[i]
          J[i, j] <- prod(sin(theta[1:(i - 1)])) * (-sin(theta[i]))
        }
      }
    }
  }

  # Last component
  for (j in seq_along(theta)) {
    if (j == 1) {
      sin_prod <- if (length(theta) > 1) {
        prod(sin(theta[2:length(theta)]))
      } else {
        1
      }
      J[n, j] <- cos(theta[1]) * sin_prod
    } else {
      prod_before <- if (j > 2) prod(sin(theta[1:(j - 1)])) else sin(theta[1])
      prod_after <- if (j < length(theta)) {
        prod(sin(theta[(j + 1):length(theta)]))
      } else {
        1
      }
      J[n, j] <- prod_before * cos(theta[j]) * prod_after
    }
  }
  J
}
