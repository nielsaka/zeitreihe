###############################################################################.
#'
#' Given some starting values, coefficients, a sequence of error vectors,
#' \code{create_varp_data} will compute a sequence of observables using a simple
#' vector autoregressive process.
#'
#' @param B A `(K x K)` matrix, providing the contemporaneous impact of the structural
#' errors on the observed variables.
#' @param W A `(K x N)` matrix, providing the sequence of structural error vectors.
#'
#' @inherit create_varp_data
#'
#'@section Details:
#'
#' * short-run restriction implemented
#'
#' @export
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

  # TODO: documentation
  # TODO: tests

  create_varp_data(A, Y0, U)
}
###############################################################################.
