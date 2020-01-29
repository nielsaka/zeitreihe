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
#' * R x K^2
#'
#' if `flatten = vec` or
#'
#' * R x (K^2 + K)/2
#'
#' if `flatten = vech`. The scalar R is the number of restrictions and K the
#' number of variables, i.e. K is the number of both the rows and the columns of
#' `X`.
#'
#' @examples
#' set.seed(8191)
#' K <- 4
#' A <- matrix(rnorm(K^2), K, K)
#' A[sample(1:K^2, K * (K-1) / 2)] <- 0
#' selection_matrix(A)
#'
selection_matrix <- function(X, value = 0:1, flatten = vec) {

  if (dim(X)[1] != dim(X)[2]) stop("X must be a square matrix.")

  X <- flatten(X)
  is_restricted <- which(X %in% value)
  R <- length(is_restricted)
  L <- length(X)

  # sometimes all elements are restricted
  # often happens in AB-model
  # if (R == L) stop("selection_matrix: all elements are restricted.")
  if (R == 0) stop("selection_matrix: no elements are restricted.")

  C <- matrix(0, R, L)
  C[cbind(seq_len(R), is_restricted)] <- 1

  return(C)
}

#' Calculate matrix ranks for identification
#'
#' There are necessary and sufficient algebraic conditions for checking whether
#' a model is (locally) identified. These conditions relate to the rank of
#' certain matrices of first-order partial derivatives. The following functions
#' compute those ranks.
#'
#' Note that `rank_AB_model` nests both `rank_A_model` and `rank_B_model` if
#' supplied with the appropriate arguments.
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
#'   involved in achieving identification. For successful identification, the
#'   rank must be equal to the number of unidentified parameters.
#'
#' * For the the A-type model, the rank has to be
#' \ifelse{latex}{\out{$K^2 + \frac{K}{2} (K + 1)$}}{K^2 + K * (K + 1) / 2}.
#'
#' * For the the B-type model, the rank has to be
#' \ifelse{latex}{\out{$K^2$}}{K^2}.
#'
#' * For the the AB-type model, the rank has to be
#' \ifelse{latex}{\out{$2K^2$}}{2K^2}.
#'
#' @examples
#' set.seed(819)
#' K <- 4
#' A <- matrix(rnorm(K^2), K, K)
#' A[sample((1:K^2)[-seq(1, K^2, K+1)], K * (K-1) / 2)] <- 0
#' diag(A) <- 1
#' A
#' SIGMA_U <- solve(A) %*% t(solve(A))
#' rank_A_model(A)
#' is_identified(A, SIGMA_U = SIGMA_U)
#'
#' K <- 3
#' A <- matrix(rnorm(K^2), K, K)
#' A[c(3, 6, 7)] <- 0
#' A[c(1, 5, 9)] <- 1
TODO: fix error message
Error in is_identified(A, SIGMA_U = SIGMA_U) : A and AB don't agree.
#'
#' @name rank_condition
NULL

#' @rdname rank_condition
rank_A_model <- function(A,
                         SIGMA_U = A_INV %*% t(A_INV),
                         C_A     = selection_matrix(A)) {
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

#' @rdname rank_condition
rank_AB_model <- function(A,
                          B,
                          SIGMA_U = A_INV %*% B %*% t(B) %*% t(A_INV),
                          C_A = selection_matrix(A),
                          C_B = selection_matrix(B)) {
  # need it? -> better in higher level function?
  # stopifnot(
  #   K == nrow(B) && K == ncol(B),
  #   K == nrow(A) && K == ncol(A),
  #   K == nrow(SIGMA_u) && isSymmetric(SIGMA_u)
  # )

  K <- var_length(A)
  D_PLUS <- duplication_matrix_ginverse(K)
  A_INV <- solve(A)

  DERIVATIVE <-
    cbind(
      rbind(
        -2 * D_PLUS %*% (SIGMA_U %x% A_INV),
        C_A,
        matrix(0, nrow = nrow(C_B), ncol = ncol(C_A))
      ),
      rbind(
        2 * D_PLUS %*% ((A_INV %*% B) %x% A_INV),
        matrix(0, nrow = nrow(C_A), ncol = ncol(C_B)),
        C_B
      )
    )

  qr(DERIVATIVE)$rank
}


# TEST
# ---

# provide A, but no SIGMA_U
# provide neither A nor B
# provide same args as for rank_AB_model
# provide non-definite / semi-definite SIGMA_U
# provide args with non matching dimensions
# make A and AB or B and AB not agree?

# TODO specify which elements are restricted? no, only 0 and 1 possible!

#' Verify whether an SVAR model is identified
#'
#' @inheritParams rank_condition
#'
#' @note
#' "The default setting assumes unit variance of the structural shocks." TRUE?
#' @return
#' @export
#'
#' @examples
is_identified <- function(A = NULL, B = NULL, SIGMA_U = NULL) {

  has_A <- !is.null(A)
  has_B <- !is.null(B)
  has_SIGMA_U <- !is.null(SIGMA_U)

  stopifnot(has_A || has_B, !has_A || has_SIGMA_U)

  if(!has_B) {
    K <- var_length(A)
    identified_A <- rank_A_model(A, SIGMA_U) == K^2 + K*(K + 1) / 2
    B <- t(chol(SIGMA_U))
  }
  if(!has_A) {
    K <- var_length(B)
    identified_B <- rank_B_model(B) == K^2
    A <- diag(K)
    SIGMA_U <- B %*% t(B)
  }

  identified <- rank_AB_model(A, B, SIGMA_U) == 2 * K^2

  if (!has_B && identified_A != identified) stop("A and AB don't agree.")
  if (!has_A && identified_B != identified) stop("B and AB don't agree.")

  return(identified)
}

##################.
#### OLD CODE ####
##################.

###############################################################################.


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
