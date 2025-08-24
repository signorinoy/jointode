# Tests for spline utility functions

test_that(".get_spline_config creates valid configuration", {
  x <- seq(0, 10, length.out = 100)

  # Test default configuration
  config <- JointODE:::.get_spline_config(x)

  expect_type(config, "list")
  expect_equal(config$degree, 3)
  expect_length(config$knots, 5)
  expect_length(config$boundary_knots, 2)
  expect_equal(config$df, 5 + 3 + 1)
  expect_equal(config$boundary_knots, range(x))

  # Test custom degree
  config2 <- JointODE:::.get_spline_config(x, degree = 2)
  expect_equal(config2$degree, 2)
  expect_equal(config2$df, 5 + 2 + 1)

  # Test custom number of knots
  config3 <- JointODE:::.get_spline_config(x, n_knots = 3)
  expect_length(config3$knots, 3)
  expect_equal(config3$df, 3 + 3 + 1)

  # Test custom boundary knots
  config4 <- JointODE:::.get_spline_config(x, boundary_knots = c(-1, 11))
  expect_equal(config4$boundary_knots, c(-1, 11))
})

test_that(".get_spline_config handles knot placement strategies", {
  x <- seq(0, 10, length.out = 100)

  # Test quantile placement (default)
  config_q <- JointODE:::.get_spline_config(
    x = x, n_knots = 4, knot_placement = "quantile"
  )
  expect_length(config_q$knots, 4)
  # Knots should be at quantiles
  expected_quantiles <- quantile(x, probs = seq(0, 1, length.out = 6)[-c(1, 6)])
  expect_equal(config_q$knots, unname(expected_quantiles))

  # Test equal spacing
  config_e <- JointODE:::.get_spline_config(
    x = x, n_knots = 4, knot_placement = "equal"
  )
  expect_length(config_e$knots, 4)
  # Check knots are equally spaced
  diffs <- diff(config_e$knots)
  expect_true(all(abs(diffs - diffs[1]) < 1e-10))

  # Test invalid placement
  expect_error(
    JointODE:::.get_spline_config(x, knot_placement = "invalid"),
    "knot_placement must be 'quantile' or 'equal'"
  )
})

test_that(".get_spline_config handles edge cases", {
  # Small data vector
  x_small <- c(1, 2, 3)
  config <- JointODE:::.get_spline_config(x_small, n_knots = 1)
  expect_length(config$knots, 1)
  expect_equal(config$knots, median(x_small))

  # Single value
  x_single <- rep(5, 10)
  config2 <- JointODE:::.get_spline_config(x_single, n_knots = 2)
  expect_length(config2$knots, 2)
  expect_true(all(config2$knots == 5))

  # With NA values
  x_na <- c(1:10, NA, NA)
  config3 <- JointODE:::.get_spline_config(x_na)
  expect_equal(config3$boundary_knots, c(1, 10))
})

test_that(".compute_spline_basis generates valid basis matrix", {
  x <- seq(0, 10, length.out = 50)
  config <- JointODE:::.get_spline_config(x, degree = 3, n_knots = 5)

  # Compute basis
  basis <- JointODE:::.compute_spline_basis(x, config)

  # Check dimensions
  expect_equal(nrow(basis), length(x))
  expect_equal(ncol(basis), config$df)

  # Check all values are finite
  expect_true(all(is.finite(basis)))

  # Check values are in reasonable range
  expect_true(all(basis >= 0))
  expect_true(all(basis <= 1))

  # Check that at least one basis function is non-zero at each point
  expect_true(all(rowSums(basis) > 0))
})

test_that(".compute_spline_basis handles extrapolation", {
  x <- seq(0, 10, length.out = 20)
  config <- JointODE:::.get_spline_config(x)

  # Points outside the boundary knots
  x_new <- c(-1, 0, 5, 10, 11)
  # Suppress expected warning about extrapolation
  basis <- suppressWarnings(JointODE:::.compute_spline_basis(x_new, config))

  # Should still return valid matrix
  expect_equal(nrow(basis), length(x_new))
  expect_equal(ncol(basis), config$df)
  expect_true(all(is.finite(basis)))
})

test_that(".compute_spline_basis_deriv generates valid derivatives", {
  x <- seq(0, 10, length.out = 50)
  config <- JointODE:::.get_spline_config(x, degree = 3, n_knots = 5)

  # Compute derivative basis
  deriv_basis <- JointODE:::.compute_spline_basis_deriv(x, config)

  # Check dimensions
  expect_equal(nrow(deriv_basis), length(x))
  expect_equal(ncol(deriv_basis), config$df)

  # Check all values are finite
  expect_true(all(is.finite(deriv_basis)))

  # Derivatives can be positive or negative
  expect_true(any(deriv_basis > 0))
  expect_true(any(deriv_basis < 0))
})

test_that(".compute_spline_basis_deriv handles degree 1 splines", {
  x <- seq(0, 10, length.out = 50)
  config <- JointODE:::.get_spline_config(x, degree = 1, n_knots = 3)

  # Compute derivative of linear splines
  deriv_basis <- JointODE:::.compute_spline_basis_deriv(x, config)

  # Should still work
  expect_equal(nrow(deriv_basis), length(x))
  expect_equal(ncol(deriv_basis), config$df)
  expect_true(all(is.finite(deriv_basis)))
})

test_that("spline basis functions have expected properties", {
  x <- seq(0, 10, length.out = 100)
  config <- JointODE:::.get_spline_config(x, degree = 3, n_knots = 5)

  basis <- JointODE:::.compute_spline_basis(x, config)

  # Each basis function should have compact support
  # (except at boundaries due to B-spline properties)
  for (j in seq_len(ncol(basis))) {
    non_zero <- which(basis[, j] > 1e-10)
    if (length(non_zero) > 0) {
      # Non-zero region should be contiguous
      expect_equal(non_zero, seq(min(non_zero), max(non_zero)))
    }
  }

  # Basis functions should sum to approximately 1 (partition of unity)
  # This is approximately true for B-splines in the interior
  interior_idx <- seq(10, 90)
  row_sums <- rowSums(basis[interior_idx, ])
  expect_true(all(abs(row_sums - 1) < 0.1))
})
