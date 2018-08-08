context("Testing SVAR data creation")

set.seed(8191)
K <- 3; N <- 1E3; p <- 2

A <- cbind(matrix(0.1, K, K), matrix(-0.05, K, K)); diag(A) <- 0.4
B <- matrix(0.4, K, K); B[upper.tri(B)] <- 0
Y0 <-matrix(0, nrow = K, ncol = p)
W <- matrix(rnorm(N * K), nrow = K, ncol = N)


test_that("Simulating SVAR works", {
  # pin down values
  out <- create_svar_data(A, B, Y0, W)
  expect_equal_to_reference(out, "svar_data.rds")

  # estimation? choleski, reduced-form
  y <- out[, (3:N) + p]
  y_1 <- out[, (3:N) + p - 1]
  y_2 <- out[, (3:N) + p - 2]

  fma <- ~ y_1[1, ] + y_1[2, ] + y_1[3, ] + y_2[1, ] + y_2[2, ] + y_2[3, ]  - 1
  resid1 <- lm(update(fma, y[1, ] ~ .))$resid
  resid2 <- lm(update(fma, y[2, ] ~ .))$resid
  resid3 <- lm(update(fma, y[3, ] ~ .))$resid

  rr <- cbind(resid1, resid2, resid3)
  covmat <- t(rr) %*% rr / (N - K)
  expect_equivalent(t(chol(covmat)), B, tol = 2E-2)
  # TODO: move to other test file and compare with chol_decomp(covmat)
})
