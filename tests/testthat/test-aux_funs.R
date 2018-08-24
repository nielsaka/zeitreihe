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


test_that("Stability criterion is check", {
  A <- matrix(c(.5, .4, .1, .5, 0, .25, 0, 0), nrow = 2)
  expect_true(check_stability(A))

  A <- matrix(c(1, 0, 0, 1), nrow = 2)
  expect_false(check_stability(A))

  A <- matrix(c(1, 0, 0, 1, -0.1, 0, 0, -0.1), nrow = 2)
  expect_true(check_stability(A))
})

test_that("Computing unconditional mean works", {
   A <- matrix(c(.5, .4, .1, .5, 0, .25, 0, 0), nrow = 2)
   nu <- c(2, 3)
   expect_equal(mean_var_process(A, nu), matrix(c(7.027027, 15.135135), 2, 1))

   A <- matrix(c(1, 0, 0, 1), nrow = 2)
   nu <- c(2, 3)
   expect_error(mean_var_process(A, nu), regexp = "process not stable")
})

test_that("Vectorising a matrix works", {
  mat <- matrix(1:8, 4, 4)
  mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]
  expect_equal(vech(mat), matrix(c(1, 2, 3, 4, 6, 7, 8, 3, 4, 8), 10, 1))
  # not identical because vech(mat) returns integer, but expected matrix with
  # 1, 2, 3, ... is numeric. would have to use 1L, 2L, ...

  mat <- matrix(1:15, 3, 5)
  vec(mat)
  expect_identical(vec(mat), matrix(1:15, 15, 1))
})


test_that("Converting to VAR(1) companion form works", {

  K <- 4
  N <- 7
  p <- 2

  A <- matrix(0.1, K, K * p)
  cf_A <- big_A(A)

  expect_identical(dim(cf_A), as.integer(c(K * p, K * p)))
  expect_identical(cf_A[1:K, ], A)
  expect_identical(
    cf_A[-(1:K), ],
    cbind(diag(K * (p - 1)), matrix(0, K * (p - 1), K))
  )

  Y <- matrix(seq_len(K * N), K, N)
  cf_Y <- big_Y(Y, p)

  expect_equal(dim(cf_Y), c(K * p, N - p + 1))
  expect_equal(c(Y[, 2], Y[, 1]), cf_Y[, 1])

  nu <- as.matrix(seq_len(K))
  cf_nu <- big_nu(nu, p)

  expect_equal(dim(cf_nu), c(K * p, 1))
  expect_equal(cf_nu[1:K, , drop = FALSE], nu)
  expect_equal(unique(cf_nu[-(1:K)]), 0)

  U <- matrix(seq_len(K * N), K, N)
  cf_U <- big_U(U, p)

  expect_equal(dim(cf_U), c(K * p, N))
  expect_equal(cf_U[1:K, ], U)
  expect_equal(unique(c(cf_U[-(1:K), ])), 0)

  SIGMA <- matrix(0.5, K, K)
  cf_SIGMA <- big_SIGMA(SIGMA, p = p)

  expect_equal(dim(cf_SIGMA), c(K * p, K * p))
  expect_equal(cf_SIGMA[1:K, 1:K], SIGMA)
  expect_equal(unique(c(cf_SIGMA[, -(1:K)], cf_SIGMA[-(1:K), ])), 0)
})

test_that("Length helpers work", {
  p <- 3
  K <- 4
  N <- as.integer(1E3)

  expect_identical(var_length(matrix(0, K, K * p)), expected = as.integer(K))
  expect_warning(var_length(1:5), "is not a matrix")

  expect_identical(lag_length(matrix(0, K, K * p)), expected = p)
  expect_error(lag_length(matrix(0, 4, 7)), "p%%1 == 0 is not TRUE")

  expect_identical(obs_length(matrix(0, K, N)), expected = N)
})
