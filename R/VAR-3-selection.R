##############################################################################.
#### Tests ----
##############################################################################.

# sequential LM tests for lag selection
VAR.Seq.test <- function(data, M, alpha) {
  K <- ncol(data)
  LR <- numeric(0)
  for (p in M:2) {
    M.large <- VAR(data, p = p)
    M.small <- VAR(data[-1, ], p = p - 1)
    LR <- c(LR, nrow(M.large$datamat) * (logLik(M.large) - logLik(M.small)))
  }
  p <- max((M:2)[LR > qchisq(1 - alpha, K ^ 2)], 1)
  list(lag.length = p, p.vals = 1 - pchisq(LR, K ^ 2))
}

# Wrap tests for serial correlation and plot p-values
SerCorr <- function(var.model, p.max, alpha = 0.05, title = "") {
  SC.pvals <- matrix(, p.max, 2)
  for (p in 1:p.max) {
    ## Portmanteau Test
    SC.pvals[p, 1] <- serial.test(var.model, lags.pt = p,
                                  type = "PT.adjusted")$serial$p.value
    ## Breusch-Godfrey Test
    SC.pvals[p, 2] <- serial.test(var.model, lags.bg = p,
                                  type = "BG")$serial$p.value
  }
  # DO NOT USE PORTMANTEAU WHEN UNCERTAIN ABOUT I(1)!!
  # plot(SC.pvals[, 1], main = "Portmanteau Test",
  # ylab = "p-values", xlab = "Horizon"); abline(h = alpha)
  plot(SC.pvals[, 2], main = title,
       ylab = "p-values", xlab = "Lag Length"); abline(h = alpha)
}
