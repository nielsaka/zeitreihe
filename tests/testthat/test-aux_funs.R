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
   nu <- as.matrix(2:3)
   expect_equal(mean_var_process(A, nu), matrix(c(7.027027, 15.135135), 2, 1))

   A <- matrix(c(1, 0, 0, 1), nrow = 2)
   nu <- as.matrix(2:3)
   expect_error(mean_var_process(A, nu), "check_stability(A) is not", f = TRUE)
})

test_that("Computing unconditional covariance works", {
  K <- 4
  p <- 2

  A <- matrix(0.0, K, K * p)
  SIGMA <- matrix(0.5, K, K)
  expect_identical(cov_var_process(A, SIGMA), SIGMA)

  p16 <- function(x) print(x, d = 16)

  A <- matrix(-0.2, K, K * p); diag(A) <- 1:K / 10
  expect_known_output(
    p16(cov_var_process(A, SIGMA)),
    "cov_var_process.txt"
  )
  expect_known_output(
    p16(cov_var_process(A, SIGMA, h = 5)),
    "cov_var_process2.txt"
  )
  expect_equal(cov_var_process(A, SIGMA, h = 150), SIGMA)
  # TODO could compare to previous code for AR(p) process
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
  set.seed(8191)

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

  expect_equal(dim(cf_Y), c(K * p, N - p))
  expect_equal(c(Y[, p + 1], Y[, p]), cf_Y[, 1])

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

  nu <- matrix(1:K, ncol = 1)
  A <- matrix(0.1, K, K * p); diag(A) <- 1:K / 10
  U <- matrix(rnorm(K * N), K, N)
  Y0 <- matrix(0, K, p)
  Y <- create_varp_data(A, Y0, U)

  cf <- companion_format(Y, nu, A, U)
  cf_U <- cf$U
  cf_U2 <- cf$Y - cf$A %*% cf$Z

  expect_equal(cf_U, cf_U2)
  expect_equal(U, selector_J(K, p) %*% cf_U2)
  expect_equal(cf$nu, matrix(c(1:K, rep(0, K)), nrow = K * p))

  # vector input fails
  expect_error(companion_format(Y[1, ], nu[1], A, U[1, ]))
  expect_error(companion_format(Y, nu, A, U[1, ]))

  # incompatible dimensions fail
  expect_error(companion_format(Y, nu, A, U[, -1]))
})

test_that("Length helpers work", {
  p <- 3
  K <- 4
  N <- as.integer(1E3)

  expect_identical(var_length(matrix(0, K, K * p)), expected = as.integer(K))
  expect_error(var_length(1:5), "is.matrix(mat) is not TRUE", f = TRUE)

  expect_identical(lag_length(matrix(0, K, K * p)), expected = p)
  expect_error(lag_length(matrix(0, 4, 7)), "p%%1 == 0 is not TRUE")
  expect_error(lag_length(1:5), "is.matrix(A) is not TRUE", f = TRUE)

  expect_identical(obs_length(matrix(0, K, N)), expected = N)
  expect_error(obs_length(1:5), "is.matrix(mat) is not TRUE", f = TRUE)
})

test_that("the sample is split", {
  K <- 1
  N <- 25
  p <- 2

  Npre <- 2
  Nest <- N
  Noos <- 5

  Nsim <- Npre + Nest + Noos

  Y <- matrix(1:(Nsim * K), K, Nsim)
  split_sample <- split_templ(Npre = Npre, Nest = Nest, Noos = Noos)
  spl <- split_sample(Y)

  mtrx <- function(x = integer(0)) matrix(x, nr = K)

  expect_type(spl, "list")
  sapply(spl, function(x) expect_is(x, "matrix"))
  expect_identical(
    object = spl,
    expected = list(
      burn_in = mtrx(),
      pre_sample = mtrx(1:2),
      estimation = mtrx(3:27),
      training   = mtrx(),
      evaluation = mtrx(),
      out_of_sample = mtrx(28:32)
    )
  )

  # wrong sample lengths: does not match ncol of Y
  split_sample_wrong <- split_templ(Npre = Npre, Nest = Nest)
  expect_error(split_sample_wrong(Y), "obs_length(Y) is not TRUE", f = TRUE)
})

test_that("Companion format is converted to original format", {
  K <- 3

  expect_is(expander_e(K), "matrix")
  expect_equal(expander_e(K, K)[1, 1], 1)
  expect_equal(unique(expander_e(K, K)[-1]), 0)
  expect_error(expander_e(K, K + 1))

  p <- 2

  expect_is(selector_J(1, 1), "matrix")
  expect_equal(selector_J(K, p)[, 1:K], diag(K))
  expect_equal(unique(c(selector_J(K, p)[, -(1:K)])), 0)
})

test_that("vech(M) is duplicated to vec(M)", {
  for (K in c(1, 10)) {
    expect_identical(dm_col_index(K), dm_col_index_2(K))

    AA <- matrix(1:K^2, K, K)
    AA[upper.tri(AA)] <- t(AA)[upper.tri(AA)]

    vec_AA <- duplication_matrix(K) %*% vech(AA)

    expect_equal(vec_AA, vec(AA))
    expect_equal(matrix(vec_AA, K, K), AA)
  }
  expect_error(duplication_matrix(5.5))
  expect_error(duplication_matrix(-1))
})
