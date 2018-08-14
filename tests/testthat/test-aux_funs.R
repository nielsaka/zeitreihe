context("Testing small auxiliary functions")

test_that("Transformation Y2Z() notation works", {

   # simple data matrix

  K <- 3
  N <- 6
  p <- 2

  Y <- matrix(seq_len(K*N), nrow = K)
  Z <- Y2Z(Y, p, const = FALSE)

  expect_equivalent(c(Y[, 2], Y[, 1]), Z[, 1])
  expect_equal(dim(Z), c(K * p, N - p))

  ###

  skip_if(save_time)

  # proper data for estimation

  K <- 2
  N <- 2E5
  p <- 3

  Y <- do.call("create_varp_data", prep_input_varp(K, N, p))
  Z <- Y2Z(Y, p, const = FALSE)
  Y <- Y[, -seq_len(p)]

  # estimate coefficients with OLS
  A <- Y %*% t(Z) %*% solve(Z %*% t(Z))

  # and compare to original values
  # deviation still quite large, but narrow down as n -> inf
  expect_equivalent(A, make_A(K, p), tolerance = 2e-3, scale = 1)
})
