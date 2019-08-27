'
# TODO
# add tests and doc for

* duplication sequence
* duplication matrix
* selection_matrix
* rank_B_model
* rank_A_model
* rank_AB_model
'

###############################################################################.
# tests
# duplication_sequence(4)

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

rank_B_model <- function(B) {
  K <- nrow(B)
  Dp <- D_plus(K)
  CB <- selection_matrix(B)
  IK <- diag(K)
  # browser()
  big <- rbind(
    2 * Dp %*% (B %x% IK),
    CB
  )
  qr(big)$rank
}

### test
# K <- 10
# all(duplication_sequence(K) == dupl_sequence_2(K))


## -----> identification !?

B = matrix(
  c(1, 0.5, 0,
    -0.2, 1, 0,
    0.8, 0, 1),
  3, 3, byrow = TRUE
)

selection_matrix <- function(Z) {
  is_zero <- vec(Z) == 0

  K <- nrow(Z)
  N <- sum(is_zero)

  res <- matrix(0, N, K^2)
  res[cbind(seq_len(N), which(is_zero))] <- 1
  return(res)
}

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



duplication_matrix_ginverse <- function(K) {
  #D <- duplication_matrix(K)
  # solve(t(D) %*% D) %*% t(D)
  MASS::ginv(duplication_matrix(K))
}

# D_plus(4)


## ----> in package "time_tools"




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


rank_B_model <- function(B) {
  K <- nrow(B)
  Dp <- D_plus(K)
  CB <- selection_matrix(B)
  IK <- diag(K)
  # browser()
  big <- rbind(
    2 * Dp %*% (B %x% IK),
    CB
  )
  qr(big)$rank
}

## -----> identification !?

B = matrix(
  c(1, 0.5, 0,
    -0.2, 1, 0,
    0.8, 0, 1),
  3, 3, byrow = TRUE
)

selection_matrix <- function(B) {
  is_zero <- vec(B) == 0

  K <- nrow(B)
  N <- sum(is_zero)

  res <- matrix(0, N, K^2)
  res[cbind(seq_len(N), which(is_zero))] <- 1
  return(res)
}
# CB <- selection_matrix(B)

D_plus <- function(K) {
  D <- duplication_matrix(K)
  solve(t(D) %*% D) %*% t(D)
}
# D_plus(4)


## ----> in package "time_tools"

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
