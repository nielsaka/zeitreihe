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
