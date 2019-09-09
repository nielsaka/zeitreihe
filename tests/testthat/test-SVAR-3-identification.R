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
  K <- 3

  ################
  # A-type model #
  ################

  required_rank_A <- K^2 + K * (K + 1) / 2

  # Cholesky
  A <- diag(K)
  A[lower.tri(A)] <- runif(K * (K - 1) / 2)

  expect_equal(rank_A_model(A), required_rank_A)

  # Non-recursive
  A <- diag(K)
  A[c(2, 6, 8)] <- runif(3)

  expect_equal(rank_A_model(A), required_rank_A)

  # insufficient number of restrictions (NOT IDENTIFIED)
  A <- diag(K)
  A[lower.tri(A)] <- runif(K * (K - 1) / 2)
  A[1, 2] <- runif(1)

  expect_lt(rank_A_model(A), required_rank_A)

  # circular (IS IDENTIFIED?)
  A <- diag(K)
  A[c(2, 6, 7)] <- runif(3)

  A_INV <- solve(A)
  SIGMA_U <- A_INV %*% diag(runif(K)) %*% t(A_INV)

  expect_lt(rank_A_model(A, SIGMA_U = SIGMA_U), required_rank_A)

  # partially circular (NOT IDENTIFIED)
  A <- diag(K)
  A[c(2, 4)] <- runif(2)

  expect_lt(rank_A_model(A), required_rank_A)

  ################
  # B-type model #
  ################

  required_rank_B <- K^2

  # Cholesky
  B <- diag(K)
  B[lower.tri(B, diag = TRUE)] <- runif(K * (K + 1) / 2)

  expect_equal(rank_B_model(B), required_rank_B)

  # Non-recursive (NOT IDENTIFIED)
  B <- diag(runif(3))
  B[c(2, 6, 8)] <- runif(3)

  expect_equal(rank_B_model(B), required_rank_B)

  # insufficient number of restrictions (NOT IDENTIFIED)
  B <- diag(K)
  B[lower.tri(B, diag = TRUE)] <- runif(K * (K + 1) / 2)
  B[1, 2] <- runif(1)

  expect_lt(rank_B_model(B), required_rank_B)

  # circular (IS IDENTIFIED?)
  B <- diag(runif(3))
  B[c(2, 6, 7)] <- runif(3)

  expect_lt(rank_B_model(B), required_rank_B)

  # partially circular (NOT IDENTIFIED)
  B <- diag(runif(3))
  B[c(2, 4)] <- runif(2)

  expect_lt(rank_B_model(B), required_rank_B)

  #################
  # AB-type model #
  #################

  required_rank_AB <- 2 * K^2

  # A Cholesky and B diagonal
  A <- diag(K)
  B <- diag(runif(K))
  A[lower.tri(A)] <- runif(K * (K - 1) / 2)

  expect_equal(rank_AB_model(A, B) - required_rank_AB,
               rank_A_model(A) - required_rank_A)

  # One restriction less (NOT IDENTIFIED)
  AA <- A
  AA[1, 3] <- runif(1)

  expect_lt(rank_AB_model(AA, B), required_rank_AB)

  BB <- B
  BB[2, 1] <- runif(1)
  expect_lt(rank_AB_model(A, BB), required_rank_AB)

  # circular (IS IDENTIFIED)

  AA[3, 1] <- 0
  expect_equal(rank_AB_model(AA, B) - required_rank_AB,
               rank_A_model(AA) - required_rank_A)

  # partially circular (NOT IDENTIFIED)

  AA[4] <- AA[6]
  AA[c(6, 7)] <- 0

  expect_lt(rank_AB_model(AA, B), required_rank_AB)

  # with one unit variance restriction, identified again

  BB <- B
  BB[1] <- 1

  expect_equal(rank_AB_model(AA, BB), required_rank_AB)

})
