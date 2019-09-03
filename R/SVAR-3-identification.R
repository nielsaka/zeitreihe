'
# add tests and doc for

* duplication sequence
* duplication matrix
* rank_B_model
* rank_AB_model
'
###########################.
#### utility functions ####
###########################.

#' Moore-Penrose inverse of duplication matrix `D`
#'
#' Simple wrapper of [MASS::ginv()].
#'
#' @param K Scalar integer, number of variables in the system.
#'
#' @return A `(K * (K + 1) / 2 x K^2)` matrix.
duplication_matrix_ginverse <- function(K) {
  MASS::ginv(duplication_matrix(K))
}

#' Selection matrix for imposing linear restrictions
#'
#' Selection matrix selecting those elements in `X` with restrictions in place.
#' This function will create the matrix `C` in `C %*% vec(X) = c`.
#'
#' @param X A square matrix, the target of imposing restrictions. The restricted
#'   elements must coincide with the values specified for `value` below. The
#'   non-restricted elements may have `NA`s or any other value except those in
#'   `value`.
#' @param value A numeric vector, the values imposed on elements of `X`. May be
#'   numeric or contain `NA`. Every element of `X` with any of these values will
#'   be assumed to be restricted to that value.
#' @param flatten A function, flattening the dimension of `X`. The default is to
#'   [vec()], but in some situations may want to select unique elements of a
#'   symmetric matrix and then [vech()] is helpful.
#'
#' @note This function is only able to create selection matrices for
#'   restrictions on individual elements. Restrictions on linear combinations of
#'   elements are not (yet) implemented.
#'
#' @return  A matrix containing 0s and 1s of dimension
#'
#' * (K^2 + K)/2 x K^2
#'
#' if `flatten = vec` or
#'
#' * (K^2 - K)/2 x (K^2 + K)/2
#'
#' if `flatten = vech`.
selection_matrix <- function(X, value = 0:1, flatten = vec) {

  if (dim(X)[1] != dim(X)[2]) stop("X must be a square matrix.")

  X <- flatten(X)
  is_restricted <- which(X %in% value)
  R <- length(is_restricted)
  L <- length(X)

  if (R == L) stop("selection_matrix: all elements are restricted.")
  if (R == 0) stop("selection_matrix: no elements are restricted.")

  C <- matrix(0, R, L)
  C[cbind(seq_len(R), is_restricted)] <- 1

  return(C)
}

#' Calculate ranks for identification
#'
#' There are necessary and sufficient algebraic conditions for checking whether
#' a model is (locally) identified. These conditions relate to the rank of
#' certain matrices of first-order partial derivatives. The following functions
#' compute those ranks.
#'
#' @param A A square numeric matrix, the coefficient of contemporaneous effects
#'   between endogenous variables. Its dimension is (K x K).
#' @param B A square numeric matrix, the coefficient of contemporaneous effects
#'   of structural shocks on endogenous variables. Its dimension is (K x
#'   K).
#' @param C_A,C_B A numeric matrix, the selection matrix for imposing linear
#'   restrictions on matrix `A` or `B`. Its dimension is (K * (K + 1) / 2 x
#'   K^2). The default assumes all elements of `A` or `B` equal to 0 or 1 are
#'   restricted to those values.
#' @param SIGMA_U A square numeric matrix, the reduced-form residual
#'   covariances. Its dimension is (K x K). The default setting assumes unit
#'   variance of the structural shocks.
#'
#' @return A scalar integer. The rank of the respective partial derivative
#'   involved in achieving identification.
#'
#' * For the the A-type model, the rank has to be
#' \ifelse{latex}{\out{$K^2 + \frac{K}{2} (K + 1)$}}{K^2 + K * (K + 1) / 2}.
#'
#' * For the the B-type model, the rank has to be
#' \ifelse{latex}{\out{$K^2$}}{K^2}.
#'
#' @name rank_condition
NULL

#' @rdname rank_condition
rank_A_model <- function(A,
                         C_A     = selection_matrix(A),
                         SIGMA_U = A_INV %*% t(A_INV)) {
  K <- var_length(A)

  D      <- duplication_matrix(K)
  D_GINV <- duplication_matrix_ginverse(K)
  A_INV  <- solve(A)
  C_s    <- selection_matrix(diag(NA, K), flatten = vech)

  DERIVATIVE <-
    cbind(
      rbind(
        - 2 * D_GINV %*% (SIGMA_U %x% A_INV),
        C_A,
        matrix(0, nrow = nrow(C_s), ncol = ncol(C_A))
      ),
      rbind(
        D_GINV %*% (A_INV %x% A_INV) %*% D,
        matrix(0, nrow = nrow(C_A), ncol = ncol(C_s)),
        C_s
      )
    )

  qr(DERIVATIVE)$rank
}

#' @rdname rank_condition
rank_B_model <- function(B, C_B = selection_matrix(B)) {

  K <- var_length(B)

  D_PLUS <- duplication_matrix_ginverse(K)
  IDENT <- diag(K)

  DERIVATIVE <- rbind(
    2 * D_PLUS %*% (B %x% IDENT),
    C_B
  )

  qr(DERIVATIVE)$rank
}

##################.
#### OLD CODE ####
##################.

###############################################################################.
# tests
# duplication_sequence(4)

###############################################################################.


is_identifiable <- function(A = NULL, B = NULL) {

  has_A <- !is.null(A)
  has_B <- !is.null(B)

  stopifnot(has_A && has_B)

  if(!has_B) {
    return(rank_A_model(A) == ...)
  }
  if(!has_A) {
    return(rank_B_model(B) == nrow(B)^2)
  }
  rank_AB_model(A, B) == ..

  # A and B
 #  check Helmut's book'
}

# NOT SURE IF AB CAN REALLY BE NESTED IN A-MODEL!
rank_AB_model <- function(A = NULL, B = NULL, SIGMA_u = NULL) {

  has_A <- !is.null(A)
  has_B <- !is.null(B)
  has_SIGMA_u <- !is.null(SIGMA_u)

  stopifnot(has_A || has_B, !has_A || has_SIGMA_u)

  if (!has_A) {
    K <- nrow(B)
    A <- diag(K)
    SIGMA_u <- diag(K) # TODO: DOES THIS MAKE SENSE?
  }
  if (!has_B) {
    K <- nrow(A)
    B <- diag(K)
  }
  if (has_A && has_B) {
    K <- nrow(A)
  }

  stopifnot(
    K == nrow(B) && K == ncol(B),
    K == nrow(A) && K == ncol(A),
    K == nrow(SIGMA_u) && isSymmetric(SIGMA_u)
  )

  D_p <- duplication_matrix_ginverse(K)

  CONTINUE:
    C_A <- selection_matrix(A) # C_A
  C_B <- selection_matrix(B)
  #C_s <- selection_matrix(SIGMA) # C_\sigma

  A_inv <- solve(A)

  BIG <- cbind(  # FIX -> use AB notation / condition
    rbind(
      -2 * D_p %*% (SIGMA_u %x% solve(A)),
      C_A,
      0
    ),
    rbind(
      Dp %*% (A_inv %x% A_inv) %*% Dp,
      0,
      C_s
    )
  )
  qr(big)$rank
}

### test
# K <- 10
# all(duplication_sequence(K) == dupl_sequence_2(K))

# CB <- selection_matrix(B)

A <- matrix(
  c(1, 0, 0,
    0, 1, 1,
    1, 1, 1),
  nrow = 3,
  ncol = 3,
  byrow = TRUE
)

# selection_matrix(A)


# D_plus(4)


##################################.
# Graphical
##################################.

# no E vertex may have same set of children

is_identifiable_graphic <- function(B) {
  # is graph disconnected?
  del <- which(rowSums(abs(B)) == 1 & colSums(abs(B)) == 1)
  if (length(del)) B <- B[-del, -del]

  K <- nrow(B)
  if (sum(B == 0) < (K * (K - 1) / 2)) return(FALSE)


  for (col in seq_len(K)) {
    i <- 1
    while(col + i <= K) {
      fail <- all((B[, col] != 0) == (B[, col + i] != 0))
      if (fail) return(FALSE)
      i <- i + 1
    }
  }
  return(TRUE)
}

# is_identifiable_graphic(B)

###############################################################################.

amat <- matrix(
  c(
    1, 0, 0.5,
    0.5, 1, 0,
    0, 0.5, 1
  ), ncol = 3, byrow = TRUE
)
# qr(amat)$rank # full rank of B matrix not sufficient!
# rank_B_model(amat)
# is_identifiable(amat)
# solve(amat)

is_identifiable <- function(B) {
  rank_B_model(B) == nrow(B)^2
}

## TEST
B = matrix(
  c(1, 0.5, 0,
    -0.2, 1, 0,
    0.8, 0, 1),
  3, 3, byrow = TRUE
)
