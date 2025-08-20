#' Validate Input Data for Joint Modeling
#'
#' @description
#' Validates longitudinal and survival data for joint ODE modeling. Checks data
#' integrity, formula validity, variable consistency, and temporal alignment
#' between longitudinal measurements and survival outcomes.
#'
#' @param longitudinal_formula A formula object specifying the longitudinal
#'   model. The left-hand side should contain the longitudinal outcome variable,
#'   and the right-hand side should specify covariates.
#' @param longitudinal_data A data frame containing longitudinal measurements
#'   with multiple observations per subject, ordered by subject ID and time.
#' @param survival_formula A formula object specifying the survival model.
#'   Must have \code{Surv()} on the left-hand side with time and
#'   status arguments.
#'   The right-hand side specifies baseline covariates.
#' @param survival_data A data frame containing survival/event data with exactly
#'   one row per subject.
#' @param id Character string specifying the name of the subject identifier
#'   variable. Must be present in both \code{longitudinal_data} and
#'   \code{survival_data}.
#' @param time Character string specifying the name of the time variable
#'   in \code{longitudinal_data} for longitudinal measurements.
#' @param spline_baseline A list of spline configuration parameters for the
#'   baseline hazard. Valid parameters include:
#'   \itemize{
#'     \item \code{degree}: Integer between 1 and 5 (default 3)
#'     \item \code{n_knots}: Integer between 1 and 20 (default 5)
#'     \item \code{knot_placement}: "quantile" or "equal" (default "quantile")
#'     \item \code{boundary_knots}: NULL or numeric vector of length 2
#'   }
#' @param spline_index A list of spline configuration parameters for the
#'   longitudinal trajectory. Valid parameters include:
#'   \itemize{
#'     \item \code{degree}: Integer between 1 and 5 (default 3)
#'     \item \code{n_knots}: Integer between 1 and 20 (default 4)
#'     \item \code{knot_placement}: "quantile" or "equal" (default "quantile")
#'     \item \code{boundary_knots}: NULL or numeric vector of length 2
#'   }
#'
#' @return Returns \code{invisible(NULL)} if validation passes.
#'   Throws an error with descriptive message if validation fails.
#'
#' @details
#' The validation process includes the following checks:
#' \itemize{
#'   \item Input type verification for all arguments
#'   \item Data dimension checks (non-empty data frames)
#'   \item Variable existence in respective data frames
#'   \item Survival formula structure (must use \code{Surv()})
#'   \item Missing value detection in critical columns
#'   \item Duplicate ID detection in survival data
#'   \item ID consistency between longitudinal and survival data
#'   \item Time value validity (non-negative for longitudinal,
#'     positive for survival)
#'   \item Temporal alignment (warns if measurements exist after
#'     event/censoring time)
#'   \item Event status validation (must be 0 or 1)
#'   \item Observation frequency assessment
#' }
#'
#' @examples
#' \dontrun{
#' # Example with simulated data
#' longitudinal_data <- data.frame(
#'   id = rep(1:100, each = 5),
#'   time = rep(0:4, 100),
#'   v = rnorm(500),
#'   x1 = rnorm(500)
#' )
#'
#' survival_data <- data.frame(
#'   id = 1:100,
#'   time = rexp(100, 0.1),
#'   status = rbinom(100, 1, 0.7),
#'   w1 = rnorm(100)
#' )
#'
#' # Validate the data
#' .validate(
#'   longitudinal_formula = v ~ x1,
#'   longitudinal_data = longitudinal_data,
#'   survival_formula = Surv(time, status) ~ w1,
#'   survival_data = survival_data,
#'   id = "id",
#'   time = "time"
#' )
#' }
#'
#' @seealso \code{\link{.process}} for data processing after validation
#'
#' @importFrom survival Surv
#' @concept utilities
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

#' Process Data for Joint Modeling
#'
#' @description
#' Transforms longitudinal and survival data into a subject-oriented structure
#' for efficient likelihood computation. Each list element contains one
#' subject's complete information including measurements, covariates, and
#' event data.
#'
#' @param longitudinal_formula A formula object specifying the longitudinal
#'   model. The response variable (left-hand side) represents the longitudinal
#'   outcome, while the right-hand side specifies time-varying covariates.
#' @param longitudinal_data A data frame containing longitudinal measurements.
#'   Should contain multiple rows per subject with measurements at
#'   different time points.
#' @param survival_formula A formula object specifying the survival model.
#'   Must use \code{Surv(time, status)} on the left-hand side.
#'   The right-hand side specifies baseline hazard covariates.
#' @param survival_data A data frame containing survival/event data.
#'   Must contain exactly one row per subject with event time and status.
#' @param id Character string specifying the name of the subject identifier
#'   variable. This variable links records between \code{longitudinal_data} and
#'   \code{survival_data}.
#' @param time Character string specifying the name of the measurement time
#'   variable in the longitudinal data.
#'
#' @return A list of class \code{"JointODE.data"} where each element
#'   represents one subject
#'   and contains:
#'   \describe{
#'     \item{id}{Subject identifier}
#'     \item{time}{Time to event or censoring}
#'     \item{status}{Event indicator (0 = censored, 1 = event)}
#'     \item{covariates}{Data frame of baseline covariates from survival model}
#'     \item{longitudinal}{List containing:
#'       \describe{
#'         \item{times}{Vector of measurement times}
#'         \item{measurements}{Vector of longitudinal outcomes}
#'         \item{covariates}{Data frame of time-varying covariates}
#'         \item{n_obs}{Number of longitudinal observations}
#'       }
#'     }
#'   }
#'   The returned list also has the following attributes:
#'   \describe{
#'     \item{n_subjects}{Total number of subjects}
#'     \item{n_observations}{Total number of longitudinal observations}
#'     \item{event_rate}{Proportion of subjects experiencing the event}
#'   }
#'
#' @details
#' This function performs the following operations:
#' \enumerate{
#'   \item Extracts model components from formulas using \code{model.frame}
#'     and \code{model.matrix}
#'   \item Reorganizes data into a subject-centric structure
#'   \item Ensures temporal ordering of longitudinal measurements
#'   \item Handles subjects with missing longitudinal data gracefully
#'   \item Computes summary statistics for the cohort
#' }
#'
#' The function is optimized for performance with pre-allocation and efficient
#' indexing strategies, making it suitable for large datasets.
#'
#' @note
#' This function assumes that \code{.validate()} has been called first to ensure
#' data integrity. Invalid data may cause unexpected behavior or errors.
#'
#' @examples
#' \dontrun{
#' # Prepare example data
#' longitudinal_data <- data.frame(
#'   id = rep(1:100, each = 5),
#'   time = rep(0:4, 100),
#'   v = rnorm(500) + rep(rnorm(100), each = 5),
#'   x1 = rnorm(500),
#'   x2 = rep(rbinom(100, 1, 0.5), each = 5)
#' )
#'
#' survival_data <- data.frame(
#'   id = 1:100,
#'   time = rexp(100, 0.1),
#'   status = rbinom(100, 1, 0.7),
#'   w1 = rnorm(100),
#'   w2 = rbinom(100, 1, 0.4)
#' )
#'
#' data_process <- .process(
#'   longitudinal_formula = v ~ x1 + x2,
#'   longitudinal_data = longitudinal_data,
#'   survival_formula = Surv(time, status) ~ w1 + w2,
#'   survival_data = survival_data,
#'   id = "id",
#'   time = "time"
#' )
#'
#' # Examine the structure
#' str(data_process[[1]]) # First subject's data
#' attributes(data_process) # Summary statistics
#' }
#'
#' @concept utilities
#'
#' @seealso
#' \code{\link{.validate}} for data validation,
#' \code{\link{JointODE}} for model fitting
#'
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

#' Configure B-Spline Parameters
#'
#' @description
#' Creates B-spline configuration with optimal knot placement for flexible
#' modeling of baseline hazards and nonlinear trajectories in joint models.
#'
#' @param x Numeric vector of data points for determining knot positions.
#'   Used to calculate quantiles when \code{knot_placement = "quantile"}.
#' @param degree Integer specifying polynomial degree of the B-spline basis.
#'   Common choices: 1 (linear), 2 (quadratic), 3 (cubic). Default is 3.
#' @param n_knots Integer specifying number of internal knots. More knots
#'   allow greater flexibility but may cause overfitting. Default is 5.
#' @param knot_placement Character string specifying knot placement strategy:
#'   \itemize{
#'     \item \code{"quantile"}: Places knots at data quantiles (default)
#'     \item \code{"equal"}: Places knots at equal intervals
#'   }
#' @param boundary_knots Numeric vector of length 2 specifying the boundary
#'   knots \code{[min, max]}, or NULL to use \code{range(x)}. Default is NULL.
#'
#' @return A list of class \code{"spline_config"} with components:
#'   \describe{
#'     \item{\code{degree}}{Polynomial degree for the B-spline basis}
#'     \item{\code{knots}}{Numeric vector of internal knot positions}
#'     \item{\code{boundary_knots}}{Numeric vector of length 2}
#'     \item{\code{df}}{Degrees of freedom (number of basis functions)}
#'   }
#'
#' @details
#' The function calculates B-spline basis parameters following the conventions
#' of \code{splines::bs()}. The degrees of freedom equals
#' \code{length(knots) + degree + 1}, which determines the number of basis
#' functions in the resulting B-spline representation.
#'
#' @note
#' Quantile-based knot placement is generally preferred as it adapts to the
#' data distribution, placing more knots where data is dense. Equal spacing
#' is useful when uniform coverage across the range is desired.
#'
#' @examples
#' \dontrun{
#' # Example 1: Configure splines for event times
#' event_times <- rexp(100, rate = 0.1)
#' config_baseline <- .get_spline_config(
#'   x = event_times,
#'   degree = 3,
#'   n_knots = 5,
#'   knot_placement = "quantile"
#' )
#'
#' # Example 2: Configure splines for index values
#' index_values <- rnorm(100)
#' config_index <- .get_spline_config(
#'   x = index_values,
#'   degree = 2,
#'   n_knots = 4,
#'   knot_placement = "equal",
#'   boundary_knots = c(-3, 3)
#' )
#' }
#'
#' @seealso \code{\link[splines]{bs}} for B-spline basis generation
#'
#' @concept utilities
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

#' Solve Joint ODE System
#'
#' @description
#' Solves the ODE system governing biomarker trajectories and cumulative hazard
#' in joint models. Computes dynamics without random effects for efficiency,
#' which are added post-computation.
#'
#' @param data A list containing subject-specific data:
#'   \describe{
#'     \item{time}{Event/censoring time}
#'     \item{status}{Event indicator (0/1)}
#'     \item{covariates}{Baseline covariates for survival model}
#'     \item{longitudinal}{List with times, measurements, covariates}
#'   }
#' @param parameters A list containing all necessary model components:
#'   \describe{
#'     \item{coef}{Model coefficients:
#'       \itemize{
#'         \item \code{baseline}: B-spline coefficients for baseline hazard
#'         \item \code{hazard}: Coefficients
#'               \eqn{[\alpha_m, \alpha_{\dot{m}}, \alpha_{\ddot{m}},}
#'               \eqn{\phi_1, ..., \phi_p]}
#'         \item \code{index_g}: B-spline coefficients for g(u) function
#'         \item \code{index_beta}: Single index coefficients
#'               \eqn{[\beta_m, \beta_{\dot{m}}, \beta_X, \beta_t]}
#'       }
#'     }
#'     \item{config}{Spline configurations:
#'       \itemize{
#'         \item \code{baseline}: Config for baseline hazard B-splines
#'         \item \code{index}: Config for single index model B-splines
#'       }
#'     }
#'   }
#' @param initial Initial conditions \eqn{[\Lambda(0), m(0), \dot{m}(0)]}.
#'   Default is c(0, 0, 0).
#'
#' @return A list containing the ODE solution with components:
#'   \describe{
#'     \item{cum_hazard}{Cumulative hazard \eqn{\Lambda(T)} at event time
#'           (without exp(b))}
#'     \item{log_hazard}{Log-hazard \eqn{\log(\lambda(T))} at event time
#'           (WITHOUT b term)}
#'     \item{biomarker}{Vector of biomarker values \eqn{m(t)} at observation
#'           times (WITHOUT b shift). NULL if no longitudinal observations
#'           exist.}
#'   }
#'
#' @details
#' The ODE system is:
#' \deqn{
#'   \frac{d}{dt}\begin{pmatrix}\Lambda \\ m \\ \dot{m}\end{pmatrix} =
#'   \begin{pmatrix}
#'     \lambda(t) \\
#'     \dot{m} \\
#'     g(\beta^T Z)
#'   \end{pmatrix}
#' }
#'
#' where:
#' \itemize{
#'   \item \eqn{\Lambda(t)} is the cumulative hazard function
#'   \item \eqn{m(t)} is the biomarker trajectory
#'   \item \eqn{\dot{m}(t)} is the biomarker velocity (first derivative)
#'   \item \eqn{\lambda(t) = \exp[\eta^T B_\lambda(t) + \alpha_m m(t) +}
#'         \eqn{\alpha_{\dot{m}} \dot{m}(t) + \alpha_{\ddot{m}} \ddot{m}(t) +}
#'         \eqn{W^T\phi]}
#'         is the hazard function
#'   \item \eqn{g(u) = \theta^T B_g(u)} is the flexible acceleration function
#'   \item \eqn{u = \beta^T Z} is the single index with
#'         \eqn{Z = [m, \dot{m}, X_{long}, t]}
#'   \item \eqn{B_\lambda(t)} and \eqn{B_g(u)} are B-spline basis functions
#' }
#'
#' The random effect \eqn{b_i} represents subject-specific deviation and is
#' excluded from the ODE for computational efficiency. To incorporate \eqn{b_i}:
#' \itemize{
#'   \item Actual log-hazard: \code{log_hazard + b}
#'   \item Actual cumulative hazard: \code{Lambda * exp(b)}
#'   \item Actual biomarker for likelihood: \code{biomarker + b}
#' }
#'
#' This separation is mathematically valid because \eqn{b_i} enters as an
#' additive constant in the log-hazard, yielding:
#' \eqn{\Lambda(t|b) = e^b \Lambda(t|0)}.
#' The ODE solver uses adaptive step size control (LSODA method) for
#' numerical stability.
#'
#' @examples
#' \dontrun{
#' # Example: Solve ODE for a subject
#' parameters <- list(
#'   coef = list(
#'     baseline = c(0.1, 0.2, 0.15), # Baseline hazard spline coefficients
#'     hazard = c(0.3, 0.1, 0.05, 0.2), # Hazard coefficients
#'     index_g = c(0.5, 0.3, 0.2), # Acceleration function coefficients
#'     index_beta = c(0.4, 0.2, 0.1, 0.15) # Single index coefficients
#'   ),
#'   config = list(
#'     baseline = baseline_spline_config,
#'     index = index_spline_config
#'   )
#' )
#'
#' # Solve the ODE system
#' sol <- .solve_joint_ode(subject_data, parameters)
#'
#' # Incorporate random effect b_i = 0.5
#' b <- 0.5
#' actual_log_hazard <- sol$log_hazard + b
#' actual_cum_hazard <- sol$cum_hazard * exp(b)
#' predicted_biomarker <- sol$biomarker + b # Add b to biomarker values
#' }
#'
#' @seealso
#' \code{\link{JointODE}} for the main model fitting function
#'
#' @concept utilities
.solve_joint_ode <- function(data, parameters, initial = c(0, 0, 0)) {
  # Helper function to compute spline basis
  .compute_spline_basis <- function(x, config) {
    splines::bs(
      x = x,
      knots = config$knots,
      degree = config$degree,
      Boundary.knots = config$boundary_knots,
      intercept = TRUE
    )
  }

  # Helper function to get longitudinal covariates at specific row
  .get_longitudinal_covariates <- function(data, row_index = NULL) {
    if (is.null(data$longitudinal$covariates) ||
          nrow(data$longitudinal$covariates) == 0) {
      return(numeric(0))
    }
    if (is.null(row_index)) row_index <- nrow(data$longitudinal$covariates)
    as.numeric(data$longitudinal$covariates[row_index, ])
  }

  # Helper function to compute acceleration
  .compute_acceleration <- function(
      biomarker, velocity, longitudinal_covariates, time_point, coef, config) {
    z <- c(biomarker, velocity, longitudinal_covariates, time_point)
    n_beta <- length(coef$index_beta)
    if (length(z) < n_beta) z <- c(z, rep(0, n_beta - length(z)))

    index_value <- sum(coef$index_beta * z[1:n_beta])
    basis_g <- .compute_spline_basis(index_value, config$index)
    as.numeric(basis_g %*% coef$index_g)
  }

  # Helper function to compute log-hazard (without random effect)
  .compute_log_hazard <- function(
      time_point, biomarker, velocity, acceleration, data, coef, config) {
    # Baseline hazard
    basis_lambda <- .compute_spline_basis(time_point, config$baseline)
    log_baseline <- as.numeric(basis_lambda %*% coef$baseline)

    # Biomarker effect
    log_biomarker <- sum(
      coef$hazard[1:3] * c(biomarker, velocity, acceleration)
    )

    # Covariate effect
    log_covariate <- 0
    if (length(coef$hazard) > 3 && !is.null(data$covariates)) {
      w <- as.numeric(data$covariates)
      if (length(w) > 0) {
        log_covariate <- sum(coef$hazard[4:(3 + length(w))] * w)
      }
    }

    log_baseline + log_biomarker + log_covariate
  }

  # Define ODE derivatives function
  ode_derivatives <- function(time, state, parameters) {
    biomarker_value <- state[2]
    biomarker_velocity <- state[3]

    # Get covariates at current time (nearest neighbor interpolation)
    has_longitudinal <- length(parameters$data$longitudinal$times) > 0
    longitudinal_covariates <- if (has_longitudinal) {
      time_index <- which.min(abs(parameters$data$longitudinal$times - time))
      .get_longitudinal_covariates(parameters$data, time_index)
    } else {
      numeric(0)
    }

    # Compute acceleration
    biomarker_acceleration <- .compute_acceleration(
      biomarker_value, biomarker_velocity, longitudinal_covariates, time,
      parameters$coef, parameters$config
    )

    # Compute hazard (without random effect b)
    log_hazard <- .compute_log_hazard(
      time, biomarker_value, biomarker_velocity, biomarker_acceleration,
      parameters$data, parameters$coef, parameters$config
    )
    hazard <- exp(log_hazard)

    # Return derivatives: [dΛ/dt, dm/dt, dṁ/dt]
    list(c(hazard, biomarker_velocity, biomarker_acceleration))
  }

  # Solve ODE System
  ode_parameters <- list(
    data = data,
    coef = parameters$coef,
    config = parameters$config
  )
  event_time <- data$time
  time_points <- sort(unique(c(0, data$longitudinal$times, event_time)))

  solution <- deSolve::ode(
    y = initial,
    times = time_points,
    func = ode_derivatives,
    parms = ode_parameters,
    method = "lsoda",
    rtol = 1e-8,
    atol = 1e-10
  )

  # Extract final state
  final_state <- solution[nrow(solution), -1]

  # Compute final quantities
  biomarker_final <- final_state[2]
  biomarker_velocity_final <- final_state[3]
  longitudinal_covariates_final <- .get_longitudinal_covariates(data)

  biomarker_acceleration_final <- .compute_acceleration(
    biomarker_final, biomarker_velocity_final, longitudinal_covariates_final,
    event_time, parameters$coef, parameters$config
  )

  # Compute log-hazard without random effect
  log_hazard_final <- .compute_log_hazard(
    event_time, biomarker_final, biomarker_velocity_final,
    biomarker_acceleration_final, data, parameters$coef, parameters$config
  )

  # Extract biomarker values at observation times
  biomarker_values <- if (length(data$longitudinal$times) > 0) {
    observation_indices <- match(data$longitudinal$times, solution[, 1])
    observed_values <- solution[observation_indices, 3] # Column 3 is m(t)
    observed_values
  } else {
    NULL
  }

  # Return results
  list(
    cum_hazard = final_state[1],
    log_hazard = log_hazard_final,
    biomarker = biomarker_values
  )
}

#' Compute Posterior Moments Using AGHQ
#'
#' @description
#' Computes posterior mean, variance, and other moments of the random effect
#' using adaptive Gauss-Hermite quadrature via the aghq package.
#'
#' @param ode_solution List containing ODE solution with components:
#'   \describe{
#'     \item{cum_hazard}{Cumulative hazard at event time without random effect}
#'     \item{log_hazard}{Log-hazard at event time without random effect}
#'     \item{biomarker}{Biomarker values at observation times without random
#'                      effect}
#'   }
#' @param data A list containing subject-specific data:
#'   \describe{
#'     \item{time}{Event/censoring time}
#'     \item{status}{Event indicator (0/1)}
#'     \item{covariates}{Baseline covariates for survival model}
#'     \item{longitudinal}{List with times, measurements, covariates}
#'   }
#' @param b_hat_init Initial guess for posterior mean \eqn{E[b_i|\mathcal{O}_i]}
#' @param measurement_error_sd Measurement error standard deviation
#' @param random_effect_sd Random effect standard deviation
#' @param k Number of quadrature points (default: 7)
#'
#' @return List containing:
#'   \describe{
#'     \item{b_hat}{Posterior mean \eqn{E[b_i|\mathcal{O}_i]}}
#'     \item{v_hat}{Posterior variance \eqn{Var[b_i|\mathcal{O}_i]}}
#'     \item{exp_b}{\eqn{E[exp(b_i)|\mathcal{O}_i]} for survival component}
#'   }
#'
#' @importFrom aghq aghq compute_moment
#'
#' @concept utilities
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

#' Compute Q Function for M-Step
#'
#' @description
#' Computes the expected complete-data log-likelihood Q(Θ|Θ^(k)) for the M-step
#' optimization in the EM algorithm. This function evaluates the objective
#' function that needs to be maximized during parameter updates.
#'
#' @param theta Vector of parameters to optimize:
#'   - baseline_spline_coefficients
#'   - hazard_coefficients
#'   - index_spline_coefficients
#'   - index_coefficients
#' @param data_list List of preprocessed subject data
#' @param posterior_list List of posterior distributions from E-step
#' @param config List containing spline configurations
#' @param sigma_e Measurement error standard deviation
#' @param sigma_b Random effect standard deviation
#' @param compute_gradient Logical, whether to compute gradient (default: FALSE)
#'
#' @return If compute_gradient = FALSE: scalar Q value
#'         If compute_gradient = TRUE: list with Q value and gradient vector
#'
#' @details
#' The Q function is decomposed as:
#' Q(Θ|Θ^(k)) = Q_long + Q_surv + Q_RE
#'
#' where:
#' - Q_long: Longitudinal component measuring fit to observations
#' - Q_surv: Survival component measuring hazard/event prediction
#' - Q_RE: Random effects component (only depends on sigma_b, not theta)
#'
#' @concept utilities
.compute_q_function <- function(
    theta, data_list, posterior_list, config,
    sigma_e, sigma_b, compute_gradient = FALSE) {

  # Parse parameter vector
  n_baseline <- config$baseline$df
  n_hazard <- length(theta) - n_baseline -
    config$index$df -
    (ncol(data_list[[1]]$longitudinal$covariates) + 3)
  n_index_g <- config$index$df

  idx <- 1
  baseline_coef <- theta[idx:(idx + n_baseline - 1)]
  idx <- idx + n_baseline
  hazard_coef <- theta[idx:(idx + n_hazard - 1)]
  idx <- idx + n_hazard
  index_g_coef <- theta[idx:(idx + n_index_g - 1)]
  idx <- idx + n_index_g
  index_beta_coef <- theta[idx:length(theta)]

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

  # Initialize gradient if needed
  if (compute_gradient) {
    grad <- numeric(length(theta))
    grad_baseline <- numeric(n_baseline)
    grad_hazard <- numeric(n_hazard)
    grad_index_g <- numeric(n_index_g)
    grad_index_beta <- numeric(length(index_beta_coef))
  }

  # Loop over subjects
  for (i in seq_along(data_list)) {
    subject_data <- data_list[[i]]
    posterior <- posterior_list[[i]]

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
        sum(residuals^2 + v_hat) / (2 * sigma_e^2) -
        subject_data$longitudinal$n_obs * log(2 * pi * sigma_e^2) / 2
      n_total_obs <- n_total_obs + subject_data$longitudinal$n_obs

      if (compute_gradient) {
        # Gradient computation for longitudinal parameters
        # This requires computing dm_i/dtheta via sensitivity ODEs
        # Placeholder: implement sensitivity analysis
      }
    }

    # Survival component
    if (subject_data$status == 1) {
      # Event occurred
      q_surv <- q_surv +
        ode_solution$log_hazard + b_hat -
        exp_b * ode_solution$cum_hazard
    } else {
      # Censored
      q_surv <- q_surv - exp_b * ode_solution$cum_hazard
    }

    if (compute_gradient) {
      # Gradient computation for survival parameters
      # This requires computing dlog_hazard/dtheta and dcum_hazard/dtheta
      # Placeholder: implement gradient computation
    }
  }

  # Random effects component (doesn't depend on theta, only on sigma_b)
  q_re <- 0
  for (i in seq_along(posterior_list)) {
    posterior <- posterior_list[[i]]
    q_re <- q_re -
      (posterior$b_hat^2 + posterior$v_hat) / (2 * sigma_b^2) -
      log(2 * pi * sigma_b^2) / 2
  }

  # Total Q function
  q_total <- q_long + q_surv + q_re

  if (compute_gradient) {
    # Combine gradients
    grad[1:n_baseline] <- grad_baseline
    grad[(n_baseline + 1):(n_baseline + n_hazard)] <- grad_hazard
    grad[(n_baseline + n_hazard + 1):
         (n_baseline + n_hazard + n_index_g)] <- grad_index_g
    grad[(n_baseline + n_hazard + n_index_g + 1):length(grad)] <-
      grad_index_beta

    return(list(value = q_total, gradient = grad))
  } else {
    return(q_total)
  }
}

#' Update Variance Components in M-Step
#'
#' @description
#' Computes closed-form updates for measurement error and random effect
#' variances in the M-step of the EM algorithm.
#'
#' @param data_list List of preprocessed subject data
#' @param posterior_list List of posterior distributions from E-step
#' @param ode_solutions List of ODE solutions for each subject
#'
#' @return List with updated variance components:
#'   \describe{
#'     \item{sigma_e}{Updated measurement error standard deviation}
#'     \item{sigma_b}{Updated random effect standard deviation}
#'   }
#'
#' @concept utilities
.update_variance_components <- function(
    data_list, posterior_list, ode_solutions) {

  n_subjects <- length(data_list)
  n_total_obs <- 0
  sum_squared_residuals <- 0
  sum_squared_random_effects <- 0

  for (i in seq_len(n_subjects)) {
    subject_data <- data_list[[i]]
    posterior <- posterior_list[[i]]
    ode_solution <- ode_solutions[[i]]

    # Measurement error variance update
    if (subject_data$longitudinal$n_obs > 0) {
      residuals <- subject_data$longitudinal$measurements -
        ode_solution$biomarker - posterior$b_hat
      sum_squared_residuals <- sum_squared_residuals +
        sum(residuals^2) +
        subject_data$longitudinal$n_obs * posterior$v_hat
      n_total_obs <- n_total_obs + subject_data$longitudinal$n_obs
    }

    # Random effect variance update
    sum_squared_random_effects <- sum_squared_random_effects +
      posterior$b_hat^2 + posterior$v_hat
  }

  # Compute updated standard deviations
  sigma_e <- sqrt(sum_squared_residuals / n_total_obs)
  sigma_b <- sqrt(sum_squared_random_effects / n_subjects)

  list(
    sigma_e = sigma_e,
    sigma_b = sigma_b
  )
}
