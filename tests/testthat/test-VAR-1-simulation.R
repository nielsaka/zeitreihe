context("Testing data creation")


test_that("AR(p) works", {
  #############################################################################.
  # basic test
  set.seed(8191)
  a <- c(0.7, -0.3, 0.2)
  ee <- rnorm(100)
  y_0 <-c(0, 0, 0)

  expect_equal_to_reference(create_arp_data(a, y_0, ee), "arp_data.rds")

  ee <- c(rnorm(100), rep(0, 1000))
  expect_equal(tail(create_arp_data(a, y_0, ee), 100), rep(0, 100))
})

test_that("VAR(1) works", {
  #############################################################################.
  # basic test
  A <- matrix(c(1/2, 1/4, 2/3, 1/3), ncol = 2)
  EE <- matrix(c(5, 4), ncol = 1)
  Y_0 <- c(10, 15)

  expect_equivalent(
    create_var1_data(A, Y_0, EE),
    matrix(c(20, 11.5), ncol = 1)
    )

  #############################################################################.
  # will it converge to uncond. mean (= zero), when stationary?
  n <- 1E4
  K <- 5
  A <- matrix(0.05, K, K)
  diag(A) <- 0.25
  EE <- matrix(0, K, n)
  EE[, seq_len(n - 1E3)] <- rnorm(K * (n - 1E3))
  Y_0 <- rep(50, K)

  expect_identical(
    create_var1_data(A, Y_0, EE)[, n], # takes nearly 1000 obs !
    c(y1 = 0, y2 = 0, y3 = 0, y4 = 0, y5 = 0)
  )
})
test_that("Var(p) works", {
  set.seed(15)
  K <- 3
  Tt <- 50
  B1 <- matrix(0.15, K, K); diag(B1) <- 0.2
  B2 <- matrix(0, K, K)
  B <- cbind(B1, B2)
  UU <- matrix(rnorm(K * Tt), K, Tt)
  Y_0 <- matrix(0, K, 1)
  Z_0 <- matrix(0, K, 2)

  out_var1 <- create_var1_data(B1, Y_0, UU[, -(1:2)])
  out_varp <- create_varp_data(B, Z_0, UU)
  expect_equal(sum(out_varp[, -(1:2)] - out_var1), 0)

  B2 <- matrix(0.05, K, K); diag(B2) <- 0.1
  B  <- cbind(B1, B2)
  out <- create_varp_data(B, Z_0, UU)
  expect_equivalent(out[1, 4, drop = FALSE], UU[, 3] %*% B[, 1] + UU[1, 4])

  #############################################################################.
  # will it converge to uncond. mean (= zero), when stationary?
  n <- 1E4
  K <- 3
  B1 <- matrix(0.2, K, K); diag(B1) <- 0.4
  B2 <- matrix(0.04, K, K); diag(B2) <- 0.09
  B <- cbind(B1, B2)
  UU <- matrix(0, K, n)
  UU[, seq_len(n - 1E3)] <- rnorm(K * (n - 1E3))
  Z_0 <- matrix(50, K, 2)

  expect_equal(
    create_varp_data(B, Z_0, UU)[, n],
    c(y1 = 0, y2 = 0, y3 = 0)
  )
})
