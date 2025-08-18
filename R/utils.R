#' Validate Input Data for Joint ODE Modeling
#'
#' @description
#' Performs comprehensive validation of longitudinal and survival data to ensure
#' they meet the requirements for joint ODE modeling. This function checks data
#' integrity, formula validity, variable consistency, and temporal alignment
#' between longitudinal measurements and survival outcomes.
#'
#' @param formula.long A formula object specifying the longitudinal model.
#'   The left-hand side should contain the longitudinal outcome variable,
#'   and the right-hand side should specify covariates.
#' @param formula.surv A formula object specifying the survival model.
#'   Must have \code{Surv()} on the left-hand side with time and
#'   status arguments.
#'   The right-hand side specifies baseline covariates.
#' @param data.long A data frame containing longitudinal measurements with
#'   multiple observations per subject, ordered by subject ID and time.
#' @param data.surv A data frame containing survival/event data with exactly
#'   one row per subject.
#' @param id Character string specifying the name of the subject identifier
#'   variable. Must be present in both \code{data.long} and \code{data.surv}.
#' @param time Character string specifying the name of the time variable
#'   in \code{data.long} for longitudinal measurements.
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
#' long_data <- data.frame(
#'   id = rep(1:100, each = 5),
#'   time = rep(0:4, 100),
#'   y = rnorm(500),
#'   x1 = rnorm(500)
#' )
#'
#' surv_data <- data.frame(
#'   id = 1:100,
#'   event_time = rexp(100, 0.1),
#'   status = rbinom(100, 1, 0.7),
#'   z1 = rnorm(100)
#' )
#'
#' # Validate the data
#' validate(
#'   formula.long = y ~ x1,
#'   formula.surv = Surv(event_time, status) ~ z1,
#'   data.long = long_data,
#'   data.surv = surv_data,
#'   id = "id",
#'   time = "time"
#' )
#' }
#'
#' @concept utilities
#'
#' @seealso \code{\link{process}} for data processing after validation
#'
#' @export
validate <- function(
    formula.long, formula.surv, data.long, data.surv, id, time) {
  # === Step 1: Validate input types ===
  stopifnot(
    "formula.long must be a formula" = inherits(formula.long, "formula"),
    "formula.surv must be a formula" = inherits(formula.surv, "formula"),
    "data.long must be a data frame" = is.data.frame(data.long),
    "data.surv must be a data frame" = is.data.frame(data.surv),
    "id must be a single character string" =
      is.character(id) && length(id) == 1,
    "time must be a single character string" =
      is.character(time) && length(time) == 1
  )

  # === Step 2: Validate data dimensions ===
  stopifnot(
    "Longitudinal data has no rows" = nrow(data.long) > 0,
    "Survival data has no rows" = nrow(data.surv) > 0
  )

  # === Step 3: Validate key variables exist ===
  # Check ID and time in longitudinal data
  if (!id %in% names(data.long)) {
    stop(sprintf("ID variable '%s' not found in longitudinal data", id))
  }
  if (!time %in% names(data.long)) {
    stop(sprintf("Time variable '%s' not found in longitudinal data", time))
  }

  # Check ID in survival data
  if (!id %in% names(data.surv)) {
    stop(sprintf("ID variable '%s' not found in survival data", id))
  }

  # === Step 4: Parse and validate survival formula ===
  surv_response <- formula.surv[[2]]
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
  long_vars <- all.vars(formula.long)
  missing_long_vars <- setdiff(long_vars, names(data.long))
  if (length(missing_long_vars) > 0) {
    stop(sprintf(
      "Variables in longitudinal formula not found in data: %s",
      paste(missing_long_vars, collapse = ", ")
    ))
  }

  # Check survival formula variables
  surv_vars <- all.vars(formula.surv)
  missing_surv_vars <- setdiff(surv_vars, names(data.surv))
  if (length(missing_surv_vars) > 0) {
    stop(sprintf(
      "Variables in survival formula not found in data: %s",
      paste(missing_surv_vars, collapse = ", ")
    ))
  }

  # === Step 6: Validate data integrity ===

  # Check for missing values in critical columns
  critical_cols <- list(
    "ID in longitudinal data" = data.long[[id]],
    "Time in longitudinal data" = data.long[[time]],
    "ID in survival data" = data.surv[[id]]
  )

  for (col_name in names(critical_cols)) {
    if (any(is.na(critical_cols[[col_name]]))) {
      stop(sprintf("Missing values found in %s", col_name))
    }
  }

  # Check for duplicate survival records
  if (any(duplicated(data.surv[[id]]))) {
    stop(paste(
      "Duplicate IDs found in survival data -",
      "each subject should have one record"
    ))
  }

  # === Step 7: Validate ID consistency ===
  long_ids <- unique(data.long[[id]])
  surv_ids <- unique(data.surv[[id]])

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
  if (any(data.long[[time]] < 0, na.rm = TRUE)) {
    stop("Negative time values found in longitudinal data")
  }

  # Check survival times are positive
  if (any(data.surv[[time_var]] <= 0, na.rm = TRUE)) {
    stop("Invalid observation times in survival data (must be positive)")
  }

  # Check for measurements beyond observation time
  subjects_with_late_obs <- character()
  for (subject_id in long_ids) {
    if (subject_id %in% surv_ids) {
      long_times <- data.long[data.long[[id]] == subject_id, time]
      surv_time <- data.surv[data.surv[[id]] == subject_id, time_var]

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
  unique_status <- unique(data.surv[[status_var]])
  invalid_status <- setdiff(unique_status, c(0, 1, NA))

  if (length(invalid_status) > 0) {
    stop(sprintf(
      "Invalid status values found: %s. Must be 0 (censored) or 1 (event)",
      paste(invalid_status, collapse = ", ")
    ))
  }

  if (any(is.na(data.surv[[status_var]]))) {
    stop("Missing values in status variable")
  }

  # === Step 10: Check observation frequency ===
  obs_per_subject <- table(data.long[[id]])
  if (all(obs_per_subject == 1)) {
    warning(paste(
      "Each subject has only one longitudinal observation -",
      "joint modeling may not be appropriate"
    ))
  }

  invisible(NULL)
}

#' Process and Reorganize Data for Joint ODE Modeling
#'
#' @description
#' Transforms validated longitudinal and survival data into an optimized list
#' structure suitable for joint ODE model fitting. Each list element contains
#' all information for a single subject, including longitudinal measurements,
#' survival outcomes, and associated covariates.
#'
#' @param formula.long A formula object specifying the longitudinal model.
#'   The response variable (left-hand side) represents the longitudinal outcome,
#'   while the right-hand side specifies time-varying covariates.
#' @param formula.surv A formula object specifying the survival model.
#'   Must use \code{Surv(time, status)} on the left-hand side.
#'   The right-hand side specifies baseline hazard covariates.
#' @param data.long A data frame containing longitudinal measurements.
#'   Should contain multiple rows per subject with measurements at
#'   different time points.
#' @param data.surv A data frame containing survival/event data.
#'   Must contain exactly one row per subject with event time and status.
#' @param id Character string specifying the name of the subject identifier
#'   variable.
#'   This variable links records between \code{data.long} and \code{data.surv}.
#' @param time Character string specifying the name of the measurement time
#'   variable
#'   in the longitudinal data.
#'
#' @return A list of class \code{"JointODE.data"} where each element
#'   represents one subject
#'   and contains:
#'   \describe{
#'     \item{id}{Subject identifier}
#'     \item{event_time}{Time to event or censoring}
#'     \item{status}{Event indicator (0 = censored, 1 = event)}
#'     \item{surv_covariates}{List of baseline covariates from survival model}
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
#' This function assumes that \code{validate()} has been called first to ensure
#' data integrity. Invalid data may cause unexpected behavior or errors.
#'
#' @examples
#' \dontrun{
#' # Prepare example data
#' long_data <- data.frame(
#'   id = rep(1:100, each = 5),
#'   time = rep(0:4, 100),
#'   y = rnorm(500) + rep(rnorm(100), each = 5),
#'   x1 = rnorm(500),
#'   x2 = rep(rbinom(100, 1, 0.5), each = 5)
#' )
#'
#' surv_data <- data.frame(
#'   id = 1:100,
#'   event_time = rexp(100, 0.1),
#'   status = rbinom(100, 1, 0.7),
#'   z1 = rnorm(100),
#'   z2 = rbinom(100, 1, 0.4)
#' )
#'
#' # First validate the data
#' validate(
#'   formula.long = y ~ x1 + x2,
#'   formula.surv = Surv(event_time, status) ~ z1 + z2,
#'   data.long = long_data,
#'   data.surv = surv_data,
#'   id = "id",
#'   time = "time"
#' )
#'
#' # Then process it
#' processed_data <- process(
#'   formula.long = y ~ x1 + x2,
#'   formula.surv = Surv(event_time, status) ~ z1 + z2,
#'   data.long = long_data,
#'   data.surv = surv_data,
#'   id = "id",
#'   time = "time"
#' )
#'
#' # Examine the structure
#' str(processed_data[[1]]) # First subject's data
#' attributes(processed_data) # Summary statistics
#' }
#'
#' @concept utilities
#'
#' @seealso
#' \code{\link{validate}} for data validation,
#' \code{\link{JointODE}} for model fitting
#'
#' @export
process <- function(
    formula.long, formula.surv, data.long, data.surv, id, time) {
  # Pre-allocate with unique IDs
  unique_ids <- unique(data.surv[[id]])
  n_subjects <- length(unique_ids)
  data.process <- vector("list", n_subjects)
  names(data.process) <- as.character(unique_ids)

  # Parse survival model once
  surv_frame <- model.frame(formula.surv, data = data.surv)
  surv_response_matrix <- model.response(surv_frame)
  has_surv_covs <- length(all.vars(formula.surv[[3]])) > 0 &&
    formula.surv[[3]] != 1
  surv_design <- if (has_surv_covs) {
    model.matrix(formula.surv, surv_frame)[, -1, drop = FALSE]
  } else {
    NULL
  }

  # Check longitudinal formula
  has_long_covs <- length(all.vars(formula.long[[3]])) > 0 &&
    formula.long[[3]] != 1

  # Build index maps for O(1) lookup
  surv_idx_map <- match(unique_ids, data.surv[[id]])
  long_id_groups <- split(seq_len(nrow(data.long)), data.long[[id]])

  for (i in seq_along(unique_ids)) {
    sid <- as.character(unique_ids[i])

    # Process longitudinal data if exists
    long_rows <- long_id_groups[[sid]]
    if (!is.null(long_rows) && length(long_rows) > 0) {
      long_subset <- data.long[long_rows, , drop = FALSE]

      # Sort once by time
      ord <- order(long_subset[[time]])
      if (!identical(ord, seq_along(ord))) {
        long_subset <- long_subset[ord, , drop = FALSE]
      }

      # Extract model components
      long_frame <- model.frame(formula.long, data = long_subset)
      long_times <- long_subset[[time]]
      long_measurements <- model.response(long_frame)
      long_covariates <- if (has_long_covs) {
        as.data.frame(
          model.matrix(formula.long, long_frame)[, -1, drop = FALSE]
        )
      } else {
        data.frame()
      }
    } else {
      long_times <- numeric(0)
      long_measurements <- numeric(0)
      long_covariates <- data.frame()
    }

    # Extract survival components (already pre-computed)
    surv_row <- surv_idx_map[i]
    event_time <- surv_response_matrix[surv_row, 1]
    event_status <- surv_response_matrix[surv_row, 2]
    surv_covariates <- if (!is.null(surv_design)) {
      as.list(surv_design[surv_row, , drop = FALSE])
    } else {
      list()
    }

    # Build subject data structure
    data.process[[i]] <- list(
      id = unique_ids[i],
      event_time = event_time,
      status = event_status,
      surv_covariates = surv_covariates,
      longitudinal = list(
        times = long_times,
        measurements = long_measurements,
        covariates = long_covariates,
        n_obs = length(long_times)
      )
    )
  }

  # Compute summary statistics efficiently
  statuses <- vapply(data.process, `[[`, numeric(1), "status")
  attr(data.process, "n_subjects") <- n_subjects
  attr(data.process, "n_observations") <- nrow(data.long)
  attr(data.process, "event_rate") <- mean(statuses == 1, na.rm = TRUE)

  data.process
}
