context("Testing SVAR identification")

test_that("utility functions work", {

  ###############################################
  # Moore-Penrose inverse of duplication matrix #
  ###############################################

  K <- 7

  D_ginv <- duplication_matrix_ginverse(K)
  D    <- duplication_matrix(K)

  expect_equal(D_ginv, solve(t(D) %*% D) %*% t(D), tol = 1E-15)
  expect_equal(ncol(D_ginv), K^2)
  expect_equal(nrow(D_ginv), (K^2 + K ) /2)

  expect_error(duplication_matrix_ginverse(0))

  ####################
  # selection matrix #
  ####################

  set.seed(8191)

  # 1. select non-restricted elements
  X <- diag(K)
  X[sample(seq_len(K^2), (K^2 - K)/2)] <- NA

  C <- selection_matrix(X)

  expect_equal(nrow(C), (K^2 + K)/2)
  expect_equal(ncol(C), K^2)

  # replace NA for computing
  X[is.na(X)] <- -999999
  expect_true(all(C %*% vec(X) >= 0))

  # 2. select elements of vech(X) below main diagonal
  X <- diag(NA, K)
  C <- selection_matrix(X, flatten = vech)

  expect_equal(nrow(C), (K^2 - K)/2)
  expect_equal(ncol(C), (K^2 + K)/2)

  # vech() checks for symmetry...
  my_vech <- function(mat) as.matrix(mat[lower.tri(mat, diag = TRUE)])

  # replace NA for computing; keep only lower triangle
  X[is.na(X) | upper.tri(X)] <- -999999
  expect_true(all(C %*% my_vech(X) == 0))

  # 3. errors
  X <- diag(K)
  expect_error(selection_matrix(X), "all elements")
  X[] <- NA
  expect_error(selection_matrix(X), "no elements")
  X <- matrix(0, 3, 2)
  expect_error(selection_matrix(X), "square matrix")
})

test_that("identification is correct", {

  set.seed(904673)

  ################
  # A-type model #
  ################

  K <- 3

  # Cholesky
  A <- diag(K)
  A[lower.tri(A)] <- runif(K * (K - 1) / 2)

  expect_equal(rank_A_model(A), K^2 + K * (K + 1) / 2)

  # Non-recursive
  A <- diag(K)
  A[c(2, 6, 8)] <- runif(3)

  expect_equal(rank_A_model(A), K^2 + K * (K + 1) / 2)

  # insufficient number of restrictions (NOT IDENTIFIED)
  A <- diag(K)
  A[lower.tri(A)] <- runif(K * (K - 1) / 2)
  A[1, 2] <- runif(1)

  expect_lt(rank_A_model(A), K^2 + K * (K + 1) / 2)

  # circular (IS IDENTIFIED?)
  A <- diag(K)
  A[c(2, 6, 7)] <- runif(3)

  A_INV <- solve(A)
  SIGMA_U <- A_INV %*% diag(runif(K)) %*% t(A_INV)

  expect_lt(rank_A_model(A, SIGMA_U = SIGMA_U), K^2 + K * (K + 1) / 2)

  # partially circular (NOT IDENTIFIED)
  A <- diag(K)
  A[c(2, 4)] <- runif(2)

  expect_lt(rank_A_model(A), K^2 + K * (K + 1) / 2)

  ################
  # B-type model #
  ################

  K <- 3

  # Cholesky
  B <- diag(K)
  B[lower.tri(B, diag = TRUE)] <- runif(K * (K + 1) / 2)

  expect_equal(rank_B_model(B), K^2)

  # Non-recursive (NOT IDENTIFIED)
  B <- diag(runif(3))
  B[c(2, 6, 8)] <- runif(3)

  expect_equal(rank_B_model(B), K^2)

  # insufficient number of restrictions (NOT IDENTIFIED)
  B <- diag(K)
  B[lower.tri(B, diag = TRUE)] <- runif(K * (K + 1) / 2)
  B[1, 2] <- runif(1)

  expect_lt(rank_B_model(B), K^2)

  # circular (IS IDENTIFIED?)
  B <- diag(runif(3))
  B[c(2, 6, 7)] <- runif(3)

  expect_lt(rank_B_model(B), K^2)

  # partially circular (NOT IDENTIFIED)
  B <- diag(runif(3))
  B[c(2, 4)] <- runif(2)

  expect_lt(rank_B_model(B), K^2)


})
