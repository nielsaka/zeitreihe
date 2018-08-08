context("Estimation routines for reduced-form VAR")

# TODO: move to test-aux_funs.R and use in other test files
# TODO: separate creation of A matrix etc. in another layer?
get_Y <- function(K, N, p) {
  set.seed(8191)
  A <- matrix(
    rep(rep(0.3^seq_len(p) * sin(pi * seq_len(p) - pi/2), each = K), times = K),
    nrow = K, ncol = K * p, byrow = TRUE)
  diag(A) <- 0.4
  Y0 <-matrix(0, nrow = K, ncol = p)
  U <- matrix(rnorm(K * N), nrow = K, ncol = N)

  create_varp_data(A, Y0, U)
}

test_that("Simple multivariate OLS succeeds", {
  p <- 2
  Y <- get_Y(K = 3, N = 1E3, p = p)
  expect_equal_to_reference(ols_mv(data = t(Y), p = p), "ols_mv.rds")

  # TODO: recover coefficients (with large data); avoid testing all the time?


})

