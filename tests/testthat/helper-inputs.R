
make_A <- function (K, p) {
  A <- rep(0.3^seq_len(p) * sin(pi * seq_len(p) - pi/2),
           each = K, times = K)
  A <- matrix(A, nrow = K, ncol = K * p, byrow = TRUE)
  diag(A) <- (0.4)
  A
}

make_Y0 <- function (K, p) {
  matrix(nu2mu(A = make_A(K, p), nu = make_nu(K)), nrow = K, ncol = p)
}

make_U <- function (K, N) {
  matrix(rnorm(K * N), nrow = K, ncol = N)
}

make_nu <- function(K) matrix(seq(5, 100, length.out = K), K, 1)

prep_input_varp <- function(K, N, p, seed = 2^13-1) {
  set.seed(seed)
  list(
    A  = make_A(K, p),
    Y0 = make_Y0(K, p),
    U  = make_U(K, N),
    nu = make_nu(K)
  )
}
