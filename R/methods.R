#' Summary of JointODE Model
#'
#' @param object A JointODE object
#' @param ... Additional arguments
#' @return A summary.JointODE object with coefficients and test statistics
#'
#' @importFrom stats coef pnorm
#' @export
summary.JointODE <- function(object, ...) {
  coefs <- coef(object)
  se <- if (!is.null(object$vcov)) {
    sqrt(diag(object$vcov))
  } else {
    rep(NA_real_, length(coefs))
  }

  structure(
    list(
      call = object$call,
      coefficients = cbind(
        Estimate = coefs,
        `Std. Error` = se,
        `z value` = coefs / se,
        `Pr(>|z|)` = 2 * pnorm(-abs(coefs / se))
      ),
      sigma = with(
        object$parameters$coefficients,
        c(sigma_e = measurement_error_sd, sigma_b = random_effect_sd)
      ),
      logLik = object$logLik,
      AIC = object$AIC,
      BIC = object$BIC,
      nobs = attr(object$data, "n_subjects"),
      convergence = object$convergence
    ),
    class = "summary.JointODE"
  )
}

#' @importFrom stats printCoefmat
#' @export
print.summary.JointODE <- function(
  x,
  digits = max(3L, getOption("digits") - 3L),
  signif.stars = getOption("show.signif.stars"),
  ...
) {
  cat("\nCall:\n")
  print(x$call)
  cat("\nVariance components:\n")
  print(x$sigma, digits = digits)
  cat("\nFixed effects:\n")
  printCoefmat(
    x$coefficients,
    digits = digits,
    signif.stars = signif.stars,
    ...
  )
  cat("\n---")
  cat("\nLog-likelihood:", x$logLik, "  AIC:", x$AIC, "  BIC:", x$BIC)
  cat("\nN =", x$nobs, " Convergence:", x$convergence$message, "\n")
  invisible(x)
}

#' Extract Model Coefficients
#'
#' @param object A JointODE object
#' @param ... Additional arguments
#' @return Named numeric vector of fixed effects coefficients
#' @export
coef.JointODE <- function(object, ...) {
  cf <- object$parameters$coefficients
  coefs <- c(cf$baseline, cf$hazard, cf$acceleration)
  names(coefs) <- c(
    paste0("baseline:", seq_along(cf$baseline)),
    paste0(
      "hazard:",
      c(
        "alpha0",
        "alpha1",
        "alpha2",
        if (length(cf$hazard) > 3) paste0("phi", seq_len(length(cf$hazard) - 3))
      )
    ),
    paste0("longitudinal:beta", seq_along(cf$acceleration))
  )
  coefs
}

#' Extract Variance-Covariance Matrix
#'
#' @param object A JointODE object
#' @param ... Additional arguments
#' @return Variance-covariance matrix of fixed effects
#' @export
vcov.JointODE <- function(object, ...) {
  object$vcov
}

#' Extract Log-Likelihood
#'
#' @param object A JointODE object
#' @param ... Additional arguments
#' @return Log-likelihood with df and nobs attributes
#' @importFrom stats coef
#' @export
logLik.JointODE <- function(object, ...) {
  structure(
    object$logLik,
    df = length(coef(object)) + 2,
    nobs = attr(object$data, "n_subjects"),
    class = "logLik"
  )
}

#' Print JointODE Model
#'
#' @param x A JointODE object
#' @param digits Number of digits for numeric output
#' @param ... Additional arguments
#' @return Invisibly returns the object
#' @importFrom stats coef
#' @export
print.JointODE <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nJoint ODE Model\n")
  cat("Call: ")
  print(x$call)
  cat(
    "\nLog-likelihood:",
    format(x$logLik, digits = digits),
    "on",
    length(coef(x)) + 2,
    "degrees of freedom\n"
  )
  cat(
    "AIC:",
    format(x$AIC, digits = digits),
    "  BIC:",
    format(x$BIC, digits = digits),
    "\n"
  )
  invisible(x)
}
