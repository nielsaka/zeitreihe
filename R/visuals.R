#' GGplot2 theme for publications
#'
#' Builds on theme_classic(). (what is modified?)
#'
#' @return A list..
#' @export
theme_publish <- function() {
  theme_classic() + #into package
    theme(
      axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5),
      axis.title = element_text(size = 16),
      text = element_text(size = 15),
      legend.key = element_rect(fill = "white", colour = "white"),
      legend.key.size = unit(1.5, "lines"),
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 16),
      legend.position = "right",
      legend.direction = "vertical",
      legend.box = "vertical",
      panel.grid.major = element_line(colour = "grey80"),
      panel.border = element_blank(),
      panel.spacing = unit(1, "lines"),
      strip.background = element_blank(),
      strip.placement = "outside"
    )
}

###############################################################################.
#' Plot Impulse Response Functions

#' @examples
#' data(oil, package = "zeitreihe")
#'
#' Y <- t(oil)
#' p <- 24
#' h <- 15
#' CI <- cholesky_irfs_wild_bootstrap(
#'   reps = 10,
#'   quantiles = c(0.16, 0.84),
#'   stdev = 1:2
#' )
#'
#' label_shocks <- c(
#'   "oil supply shock",
#'   "aggregate demand shock",
#'   "oil specific-demand shock"
#' )
#'
#' # normalise shocks such that response of first variable is negative on impact
#' norm <- matrix(c(0, 0, 1, 0, 0, 1, 0, 0, 1), 3, 3)
#'
#' # cumulate the response of the first variable
#' cumulate <- 1
#'
#' IRF <- cholesky_irfs(Y, p, h, CI)
#' IRF[IRF$shock == "oil supply shock", ]$value <- - IRF[IRF$shock == "oil supply shock", ]$value
#'
#'plot_irfs(IRF)

# TODO: print legend?
# TODO: document -> especially 'lty'

plot_irfs <- function(IRF, lty) {

  # default lines
  if (missing(lty)) {
    lty <- rbind(
      c("sdp1", "dashed"),
      c("sdm1", "dashed"),
      c("sdp2", "longdash"),
      c("sdm2", "longdash"),
      c("q0.025", "twodash"),
      c("q0.05", "dotdash"),
      c("q0.16", "dotted"),
      c("q0.84", "dotted"),
      c("q0.95", "dotdash"),
      c("q0.975", "twodash")
    )
    colnames(lty) <- c("stat", "linetype")
  }

  ci_stats <- unique(IRF$stat)
  ci_stats <- ci_stats[ci_stats != "point"]

  if (!all(ci_stats %in% lty[, "stat"])) {
    stop("Unknown stats. Configure custom statistic and linetype via 'lty'.")
  }

  if (!require("ggplot2")) {
    stop("Package \"pkg\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  irf_vis <-
    ggplot(IRF) +
    geom_line(aes(x = h, y = value), data = IRF[IRF$stat == "point", ]) +
    geom_hline(yintercept = 0) +
    facet_grid(response ~ shock, scales = "free", switch = "y") +
    theme_publish() +
    ylab("") +
    xlab("Horizon")

  for (stat in ci_stats) {
    irf_vis <-
      irf_vis +
      geom_line(
        aes(x = h, y = value),
        data = IRF[IRF$stat == stat, ],
        linetype = lty[lty[, "stat"] == stat, "linetype"]
      )
  }
  plot(irf_vis)
}
