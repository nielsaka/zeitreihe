

#' A <- matrix(c(0, 0, 1, 0, 0, 1, 0, 1, 0), 3, 3)
find_ident_node <- function(A) {

  A & t(A)

  ident_node <- rowSums(A) == 0
}

find_instruments <- function(adjmat) {

  K <- nrow(A)
  A <- adjmat

  ident_nodes <- rowSums(A) == 0
  #ident_nodes <- which(rowSums(A) == 0)

  # use possible instruments
  AArem <- AA <- A
  AArem[, ident_nodes] <- 0

  # achieves identification?
  useful <- which(rowSums(AArem & t(AArem))[ident_nodes] != 0)
  if (length(useful)) {
    # record instruments; could be several..

  } else {
    # check one edge further
    AArem <- AA <- AA %*% A


    AAmir <- AA & t(AA)

    #
    AA[, ident_nodes, drop = FALSE]

    #

    ident_nodes %*% A

    A[, ident_nodes] %*% A

  }



  find_ident_node(A)


}



# restriction on A possible -> generalise by function that accepts input?
rand_coeff_iv <- function(
  random,
  draw_theta = runif(L, -1, 1),
  restriction = NULL
) {
  stopifnot(nrow(random) == ncol(random))

  L <- sum(random)
  K <- ncol(random)


  is_restricted <-
    if (is.null(restricted)) matrix(0, K, K) else !is.na(restriction)

  src_iv_adj_mat <- matrix(1, K, K) - diag(K) - random - is_restricted

  find_instruments(src_iv_adj_mat)


}

IV <- function(y, X, Z) {
  y %*% t(Z) %*% solve(X %*% t(Z))
 # solve(t(Z) %*% X) %*% t(Z) %*% y
}

########################################################
### Hard coded example from Kilian and Murphy (2012) ###
########################################################

#' @examples
#' data(oil, package = "zeitreihe")
#' rf_model <- ols_mv(t(oil), p = 24)
#' (B <- rand_coeff_kilian(rf_model))
rand_coeff_kilian <- function(rf_model) {
  Y <- rf_model$Y
  p <- lag_length(rf_model$BETA.hat[, -1])
  Z <- Y2Z(Y, p)
  YN <- Y[, -seq_len(p), drop = FALSE]

  # NOTE: SIGN FLIPS BELOW AGAIN!?
  theta <- runif(3, -1, 1)
  B <- diag(3)
  B[cbind(c(1, 1, 2), c(2, 3, 1))] <- theta / (1 - abs(theta))

  eps_bar1 <- B[1, ] %*% YN

  # y is regressand
  # X are regressors
  # Z are instruments

  # order of instruments does not matter!?

  # second equation
  yN <- B[2, ] %*% YN
  XN <- rbind(YN[3, ], Z)
  ZN <- rbind(eps_bar1, Z)

  B[2, 3] <- - IV(yN, XN, ZN)[1]
  eps_bar2 <- B[2, ] %*% YN

  # third equation
  yN <- B[3, ] %*% YN
  XN <- rbind(YN[1:2, ], Z)
  ZN <- rbind(eps_bar1, eps_bar2, Z)

  B[3, 1:2] <- - IV(yN, XN, ZN)[1:2]

  # NO!
  # K <- var_length(B)
  # B[!diag(K)] <- -1 * B[!diag(K)]

  solve(B) # effect of one standard deviation shocks
}

#' @examples
#' data(oil, package = "zeitreihe")
#' rf_model <- ols_mv(t(oil), p = 24)
#'
#' set.seed(1345)
#' (B <- rand_coeff_kilian(rf_model))
draw_sign_kilian <- function(rf_model) {

  signs <- matrix(c(-1, -1, 1, 1, 1, 1, 1, -1, 1), 3, 3)

  reps <- 0
  sign_r <- 1

  while (!all(sign_r == 0 | sign_r == 3)) {
    B <- rand_coeff_kilian(rf_model)
    sign_r <- colSums(B * signs > 0)
    reps <- reps + 1
  }
  # may need to flip signs of column
  flip <- sign_r == 0
  B[, flip] <- -1 * B[, flip]

  # names
  dimnames(B) <- list(
    c("production", "activity", "price"),
    c("supply", "aggregate demand", "oil-specific demand")
  )

  list(B = B, reps = reps)
}

#' @examples
#' data(oil, package = "zeitreihe")
#' rf_model <- ols_mv(t(oil), p = 24)
#'
#' set.seed(1345)
#' B_inv <- ouliaris_pagan_sign(rf_model, total_models = 10)
ouliaris_pagan_sign <- function(rf_model, total_models) {
                                #rand_coeff, check_sign, reps) {
  reps <- models <-  vector("list", total_models)
  for(i in seq_len(total_models)) {
    out <- draw_sign_kilian(rf_model)
    models[[i]] <- out$B
    reps[[i]] <- out$reps
    cat("Draw: model number", i, "out of", total_models, "\n")
  }
  list(B = models, reps = reps)
}

#' @examples
#' data(oil, package = "zeitreihe")
#' rf_model <- ols_mv(t(oil), p = 24)
#'
#' set.seed(1345)
#' BEES <- ouliaris_pagan_sign(rf_model, total_models = 10)
#' IRF_sign <- IRF_multimod(rf_model, BEES)
IRF_multimod <- function(rf_model, BEES, h = 18) {
  PHI <- MA_coeffs(rf_model$BETA.hat[, -1], h)

  cols <- c("response", "shock", "h", "value", "model")
  IRF <- data.frame(matrix(nc = length(cols), nr = 0, dim = list(NULL, cols)))

  for (i in seq_along(BEES$B)) {
    C <- sMA_coeffs(PHI, BEES$B[[i]])
    # cumulate effect on production
    C[1, , ] <-
      aperm(apply(C[1, , , drop = FALSE], c(1, 2), cumsum), c(2, 3, 1))
    C <- as.data.frame.table(C)
    C$model <- i
    colnames(C) <- c("response", "shock", "h", "value", "model")
    IRF <- rbind(IRF, C)
  }
  IRF
}

if (FALSE) {
  require(ggplot2)

  data(oil, package = "zeitreihe")
  rf_model <- ols_mv(t(oil), p = 24)

  # set.seed(1345)
  set.seed(145)
  tictoc::tic()
  BEES <- ouliaris_pagan_sign(rf_model, total_models = 500)
  tictoc::toc()
  IRF_sign <- IRF_multimod(rf_model, BEES)


  irf_vis <-
    ggplot(IRF_sign) +
    geom_line(aes(x = h, y = value, group = model)) +
    geom_hline(yintercept = 0) +
    facet_grid(response ~ shock, scales = "free", switch = "y") +
    theme_publish() +
    ylab("") +
    xlab("Horizon")

plot(irf_vis)
}

###########################################################
### Hard code example with additional zero restrictions ###
###########################################################

#' @examples
#' data(oil, package = "zeitreihe")
#' rf_model <- ols_mv(t(oil), p = 24)
#' (B <- rand_coeff_kilian(rf_model))
rand_coeff_zero_kilian <- function(rf_model) {
  Y <- rf_model$Y
  p <- lag_length(rf_model$BETA.hat[, -1])
  Z <- Y2Z(Y, p)
  YN <- Y[, -seq_len(p), drop = FALSE]

  # NOTE: SIGN FLIPS BELOW AGAIN!?
  theta <- runif(1, -1, 1)
  B <- diag(3)
  B[1, 3] <- theta / (1 - abs(theta))

  eps_bar1 <- B[1, ] %*% YN

  # y is regressand
  # X are regressors
  # Z are instruments

  # order of instruments does not matter!?

  # second equation
  yN <- B[2, ] %*% YN
  XN <- rbind(YN[3, ], Z)
  ZN <- rbind(eps_bar1, Z)

  B[2, 3] <- - IV(yN, XN, ZN)[1]
  eps_bar2 <- B[2, ] %*% YN

  # third equation
  yN <- B[3, ] %*% YN
  XN <- rbind(YN[1:2, ], Z)
  ZN <- rbind(eps_bar1, eps_bar2, Z)

  B[3, 1:2] <- - IV(yN, XN, ZN)[1:2]

  # NO!
  #K <- var_length(B)
  #B[!diag(K)] <- -1 * B[!diag(K)]

  solve(B) # effect of one standard deviation shocks
}

#' @examples
#' data(oil, package = "zeitreihe")
#' rf_model <- ols_mv(t(oil), p = 24)
#'
#' set.seed(1345)
#' B_inv <- ouliaris_pagan_sign(rf_model, total_models = 10)
ouliaris_pagan_sign_zero <- function(rf_model, total_models) {
  #rand_coeff, check_sign, reps) {
  reps <- models <-  vector("list", total_models)
  for(i in seq_len(total_models)) {
    out <- draw_sign_zero_kilian(rf_model)
    models[[i]] <- out$B
    reps[[i]] <- out$reps
    cat("Draw: model number", i, "out of", total_models, "\n")
  }
  list(B = models, reps = reps)
}

#' @examples
#' data(oil, package = "zeitreihe")
#' rf_model <- ols_mv(t(oil), p = 24)
#'
#' set.seed(1345)
#' (B <- rand_coeff_kilian(rf_model))
draw_sign_zero_kilian <- function(rf_model) {

    signs <- matrix(c(-1, -1, 1, 1, 1, 1, 1, -1, 1), 3, 3)

    reps <- 0
    sign_r <- 1

    while (!all(sign_r == 0 | sign_r == 3)) {
      B <- rand_coeff_zero_kilian(rf_model)
      sign_r <- colSums(B * signs > 0)
      reps <- reps + 1
    }
    # may need to flip signs of column
    flip <- sign_r == 0
    B[, flip] <- -1 * B[, flip]

    # names
    dimnames(B) <- list(
      c("production", "activity", "price"),
      c("supply", "aggregate demand", "oil-specific demand")
    )

    list(B = B, reps = reps)
  }
