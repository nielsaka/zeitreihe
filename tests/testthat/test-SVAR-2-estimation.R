context("Testing SVAR estimation")

set.seed(8191)

# number of variables, observations and lag length
K <- 3
N <- 1E6
p <- 2

# prepare input
A <- cbind(matrix(0.1, K, K), matrix(-0.05, K, K)); diag(A) <- 0.4
Be <- matrix(0.4, K, K); Be[upper.tri(Be)] <- 0; diag(Be) <- 1
Y0 <-matrix(0, nrow = K, ncol = p)
W <- matrix(rnorm(N * K), nrow = K, ncol = N)

# draw data
Y <- create_svar_data(A, Be, Y0, W)

# estimate with OLS
Be_hat_ols <- ols_cholesky(Y, p)

test_that("Covariance matrix is decomposed via Cholesky", {
  expect_equivalent(Be_hat_ols, Be, tol = 2E-3)
})

# set restrictions
By_init <- diag(K)
Be_init <- matrix(0, K, K)
Be_init[lower.tri(Be_init, diag = TRUE)] <- NA

test_that("concentrated log-lik fn is initialised and evaluated", {
  log_lik <- conc_log_lik_init(Y, p, By_init, Be_init)
  expect_equal(log_lik(rep(0.35, 6)), -14793189.4491467550396919)
})

# estimate with MLE
Be_hat_mle <- mle_svar(Y, p, By_init, Be_init)$Be

test_that("concentrated log-lik of SVAR is maximised", {
  expect_known_output(print(Be_hat_mle, d = 16), "Be_hat_mle.txt")
  expect_equivalent(Be_hat_mle, Be, tol = 2E-3)
  # TODO how precise should MLE be?? Can it be more precise? Gradient? Optimiser?
  expect_equivalent(Be_hat_mle, Be_hat_ols, tol = 1E-6)
})

##################.
skip_if(save_time)
##################.

# valid over-identifying restriction
By <- diag(K)
Be_init[3, 3] = 1

# unrestricted
ols_fit <- ols_mv(Y, p)
mle_fit <- mle_var(Y, p) ### takes ages .... !

# restricted
mle_fit_res <- mle_svar(Y, p, By = By_init, Be = Be_init)
Be_hat_mle_res <- mle_fit_res$Be
By_hat_mle_res <- mle_fit_res$By
By_hat_mle_res_inv <- solve(By_hat_mle_res)
SIGMA_r <- By_hat_mle_res_inv %*% Be_hat_mle_res %*% t(Be_hat_mle_res) %*% t(By_hat_mle_res_inv)

# unrestriced
SIGMA_mle_ur <- mle_fit$SIGMA.hat
SIGMA_ols_ur <- ols_fit$SIGMA.hat

# TODO expectation
SIGMA_ols_ur - SIGMA_mle_ur

# rescale OLS estimate to obtain ML estimate
Be_hat_mle_precise <- chol_decomp(ols_fit$SIGMA.hat * (N - (K * p) - 1) / N)

# TODO expecation
Be_hat_mle_res - Be_hat_mle_precise

SIGMA_true <- tcrossprod(Be)

# TODO expecation
SIGMA_ols_ur - SIGMA_true
SIGMA_mle_ur - SIGMA_true
SIGMA_r - SIGMA_true

# TODO expecation
SIGMA_r - SIGMA_ols_ur
SIGMA_r - SIGMA_mle_ur

# TODO expecation
log(det(SIGMA_r)) - log(det(SIGMA_mle_ur))

z_r <- determinant(SIGMA_r)
ldSIGMA_r <- as.numeric(z_r$sign * z_r$modulus)

z_ur <- determinant(SIGMA_mle_ur)
ldSIGMA_mle_ur <- as.numeric(z_ur$sign * z_ur$modulus)

# TODO expecation
## Likelihood ratio test
N * (ldSIGMA_r - ldSIGMA_mle_ur)

# TODO
# numerical precision; bad optimisation? better with gradient provided?
# should really be very close to 0
