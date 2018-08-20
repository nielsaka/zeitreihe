###############################################################################.
#' Create data using a structural VAR
#'
#' Given some starting values, coefficients, a sequence of *structural* error
#' vectors, and a matrix of contemporaneous impact effects of those errors,
#' \code{create_varp_data} will compute a sequence of observables using a simple
#' vector autoregressive process.
#'
#' @param B A `(K x K)` matrix, providing the contemporaneous impact of the
#'   structural errors on the observed variables.
#' @param W A `(K x N)` matrix, providing the sequence of structural error
#'   vectors.
#'
#' @inherit create_varp_data
#' @family functions for creating data
#'
#' @section Details:
#' * short-run restriction implemented
#' * could just pass structural errors to creat_varp_data; however, here, it's
#' documented.
#'
#' @export
#' @examples
#' K <- 3
#' N <- 1E3
#' p <- 2
#'
#' A <- cbind(matrix(0.1, K, K), matrix(-0.05, K, K)); diag(A) <- 0.4
#' B <- matrix(0.4, K, K); B[upper.tri(B)] <- 0
#' Y0 <-matrix(0, nrow = K, ncol = p)
#'set.seed(8191)
#' W <- matrix(rnorm(N * K), nrow = K, ncol = N)
#'
# Y <- create_svar_data(A, B, Y0, W)
# TODO no constant yet ?!
create_svar_data <- function(A, B, Y0, W) {
  # B - contemporaneous impact of structural innovations
  # A - reduced-form lag coefficients
  # C - structural lag coefficients
  # U - reduced-form innovations
  # W - structural innovations

  # turn into reduced form
  U <- B %*% W

  # normalise errors?
  # U <- U / sqrt(diag(B %*% t(B)))

  # TODO: tests

  create_varp_data(A, Y0, U)
}
###############################################################################.
