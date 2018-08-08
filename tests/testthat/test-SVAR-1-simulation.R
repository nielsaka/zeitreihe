context("Testing SVAR data creation")

set.seed(8191)
K <- 3; N <- 1E3; p <- 2

A <- cbind(matrix(0.1, K, K), matrix(-0.05, K, K)); diag(A) <- 0.4
B <- matrix(0.4, K, K); B[upper.tri(B)] <- 0
Y0 <-matrix(0, nrow = K, ncol = p)
W <- matrix(rnorm(N * K), nrow = K, ncol = N)


test_that("Simulating SVAR works", {
  # pin down values
  expect_equal_to_reference(create_svar_data(A, B, Y0, W), "svar_data.rds")

  # estimation? choleski, reduced-form


})
