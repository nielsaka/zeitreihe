###############################################################################.
#'Transform a matrix of time-series by adding lags
#'
#'A matrix Y is transformed into Z notation (see Lütkepohl (2005, p. 70)). This
#'facilitates estimation of VAR(p) models.
#'
#'@param Y A (K x T) matrix, containing the data of the K variables.
#'@param p An integer, the number of lags to create. Must be greater or equal 0.
#'
#'@return A ((K * (p + 1)) x (T - p)) matrix, the first K rows will contain the
#'  original data and the remaining rows will contain its lags up to order p.
#'  The columns will contain time periods, of which there are now T - p left
#'  after setting aside pre-sample values for the lags.
#'
#'@section Note: For regression purposes, the first K rows would have to be
#'  discarded as they contain Y itself; see example.
#'
#'@examples
# TODO: example not working anymore!
#' Y <- draw_data_DH08_m1(n = 5e5)
#'
#' X <- Y2Z(Y, p = 4)
#' Y <- X[seq_len(4), ] # has correct sample size
#' Z <- X[-seq_len(4), ]
#'
#' B <- Y %*% t(Z) %*% solve(Z %*% t(Z))
#'
#'@section ToDo:
#' add sanity checks... data.frame? matrix? dimension? max, min p?
#' no intercept (yet)!
#' rownames?
#'
Y2Z <- function(Y, p) {
  # TODO: add sanity checks... data.frame? matrix? dimension?
  # no intercept (yet)!
  Y <- t(as.matrix(Y))
  t(embed(Y, p + 1))
}
