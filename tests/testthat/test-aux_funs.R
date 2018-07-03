context("Testing small auxiliary functions")

test_that("Transformation Y2Z() notation works", {
  # simple data matrix with K = 5 variables and T = 6 observations, p = 3 lags.
  K <- 5; T <- 6; p <- 3
  Y <- t(matrix(seq_len(K*T), ncol = K))
  out <- Y2Z(Y, p)

  expect_equivalent(Y[1, -seq_len(p)], out[1, ])
  expect_equivalent(Y[3, seq_len(p)], out[K*p + 3, ])

  # proper data and use for estimation
  set.seed(2^12-1)
  K <- 3; Tt <- 5E4
  B1 <- matrix(0.15, K, K); diag(B1) <- 0.2
  B2 <- matrix(0, K, K)
  B <- cbind(B1, B2)
  UU <- matrix(rnorm(K * Tt), K, Tt)
  Z_0 <- matrix(0, K, 2)

  Y <- create_varp_data(B, Z_0, UU)

  X <- Y2Z(Y, p = 4)
  Y <- X[seq_len(K), ]
  Z <- X[-seq_len(K), ]

  # Estimate coefficients with OLS
  B3 <- Y %*% t(Z) %*% solve(Z %*% t(Z))

  # and compare to original values
  # deviation still quite large, but narrow down as n -> inf
  expect_equal(sum(abs(B3)), sum(abs(B)), tolerance = 15e-2, scale = 1)
})
