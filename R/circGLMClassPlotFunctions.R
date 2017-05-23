#' Plot circGLM Object
#'
#' General plot function for circGLM objects, which dispatches the chosen type
#' of plotting to the corresponding function.
#'
#' @param m A circGLM object to be plotted.
#' @param type Character string giving the type of plotting. The options are
#'   \code{"trace"}, \code{"tracestack"}, \code{"predict"}, \code{"meancompare"}
#'   and \code{"meanboxplot"}.
#' @param ... Additional arguments to be passed to print functions.
#'
#' @export
#'
#' @seealso \code{\link{"plot_trace.circGLM"}},
#'   \code{\link{"plot_tracestack.circGLM"}},
#'   \code{\link{"plot_predict.circGLM"}},
#'   \code{\link{"plot_meancompare.circGLM"}} and
#'   \code{\link{"plot_meanboxplot.circGLM"}}.
#'
#' @examples
#' plot(circGLM(rvmc(10, 1, 1)))
#'
#' dat <- generateCircGLMData(n = 100, nconpred = 1, ncatpred = 1)
#' m   <- circGLM(th = dat[, 1], X = dat[, -1])
#'
#' # Traceplot by default
#' plot(m)
#'
#' # Traceplot stack
#' plot(m, type = "tracestack")
#'
#' # Prediction plot
#' plot(m, type = "predict")
#'
#' # Mean comparisons
#' plot(m, type = "meancompare")
#' plot(m, type = "meanboxplot")
#'
plot.circGLM <- function(m, type = "trace", ...) {
  printFunName <- paste0("plot_", type, ".circGLM")
  do.call(printFunName, args = c(list(m = m), list(...)))
}



#' Create a prediction plot from circGLM objects
#'
#' Plot the predictions made by a circGLM analysis.
#'
#' Creates a ggplot showing a prediction plot showing linear predictor against
#' the circular outcome, an with an optional grouping variable. One or more
#' regression lines show the predicted values for different values of the linear
#' and categorical predictors.
#'
#' Predictors \code{x} and \code{d} and outcome \code{th} can be provided as
#' numeric vectors of the same length as the outcome in the circGLM object
#' \code{m}. This allows plotting the regression line from an earlier dataset on
#' a new dataset.
#'
#' Alternatively, \code{x} and \code{d} can be strings containing names of
#' corresponding predictors in the original model. In that case, \code{th}
#' should not be provided.
#'
#' The function makes an effort to find predictors to plot if none are given,
#' where it will simply take the first predictor in the dataset. If a plot
#' without grouping is required, \code{d} can be set to \code{NA}.
#'
#'
#' @param m A circGLM object.
#' @param x Optional; Either a numeric vector with a continuous predictor or
#'   string naming the desired variable to plot on the x-axis. If missing, we
#'   just use the first continuous predictor in the circGLM object.
#' @param d Optional; Either a numeric vector with a categorical predictor or
#'   string naming the desired variable to plot on the x-axis. If missing, we
#'   just use the first categorical predictor in the circGLM object.
#' @param th Optional; Can be a new numeric vector containing outcome angles
#'   corresponding to predictors \code{x} and potentially \code{d}.
#' @param linkfun The link function to be used. Should be the same as was used
#'   for the creation of the circGLM object.
#' @param xlab A character string with the x-label.
#' @param ylab A character string with the y-label.
#' @param colorPalette The colors to use in plotting, max 2.
#'
#' @return A \code{\link[ggplot2]{ggplot2}} plot, to which further \code{ggplot}
#'   elements can be added.
#' @export
#'
#' @seealso \code{\link{"plot_trace.circGLM"}},
#'   \code{\link{"plot_tracestack.circGLM"}},
#'   \code{\link{"plot_meancompare.circGLM"}},
#'   \code{\link{"plot_meanboxplot.circGLM"}},
#'   \code{\link{"plot.circGLM"}}.
#'
#' @examples
#' plot(circGLM(rvmc(10, 1, 1)))
#'
#' dat <- generateCircGLMData()
#' plot(circGLM(th = dat[, 1], X = dat[, -1]), type = "predict")
#'
plot_predict.circGLM <- function(m, x, d, th,
                                 linkfun = function(x) atanLF(x, 2),
                                 xlab = NA, ylab = expression(theta),
                                 colorPalette = c("#E69F00", "#56B4E9")) {

  # Try to get a label for x.
  if (is.na(xlab) && is.character(x)) xlab <- x else xlab <- "x"

  # Usually, we have no th given, so we use the outcome from the model m.
  if (missing(th)) th <- m$data_th

  # First, try to find x and d where possible.
  if (missing(x)) {
    # If we don't have a given x, just use the first one if there is one.
    if (ncol(m$data_stX) != 0) {
      btx <- m$bt_mean[,1]
      x   <- m$data_stX[, 1, drop = FALSE]
    } else {
      warning(paste("Predict plot can not be drawn without continuous",
                    "predictors. Returning a mean posterior plot instead."), call. = FALSE)
      return(plot_meancompare.circGLM(m))
    }
    # Find the correct column x if it is given as a string.
  } else if (is.character(x)) {
    btx <- m$bt_mean[, x]
    x   <- m$data_stX[, x, drop = FALSE]
  } else {
    btx <- m$bt_mean[, colnames(x)[1]]
  }

  if (missing(d)) {
    # If we don't have a given d, just use the first one if there one.
    if (ncol(m$data_d) != 0) {
      d   <- m$data_d[, 1, drop = FALSE]
      dtd <- m$dt_meandir[, 1]
      pdat <- data.frame(th = th, x = x, d = d)
    } else {
      # In this case there is no grouping.
      pdat <- data.frame(th = th, x =  x)
      d <- NA
    }
    # Find the correct column d if it is given as a string.
  } else if (is.character(d)) {
    d   <- m$data_d[, d]
    dtd <- m$dt_meandir[, d]
    pdat <- data.frame(th = th, x =  x, d = d)
  } else if (is.na(d)) {
    # In this case there is no grouping.
    pdat <- data.frame(th = th, x =  x)
  } else {
    # In the final alternative d has been provided.
    dtd <- m$dt_meandir[, colnames(d)[1]]
  }

  # The base prediction function.
  predfun <- function(x) m$b0_meandir + linkfun(x * btx)



  # Check if there is a grouping, then return the appropriate plot.
  if ((is.na(d) | missing(d)) || ncol(d) == 0) {
    p <- ggplot2::ggplot(data = pdat, ggplot2::aes(y = th, x = x)) +
      ggplot2::geom_point() +
      ggplot2::stat_function(fun = predfun,
                             col = colorPalette[1]) +
      ggplot2::theme_bw() +
      ggplot2::scale_colour_manual(values = colorPalette) +
      ggplot2::ylab(ylab) +
      ggplot2::xlab(xlab)
  } else {
    p <- ggplot2::ggplot(data = pdat, ggplot2::aes(y = th, x = x, col = factor(d))) +
      ggplot2::geom_point() +
      ggplot2::stat_function(fun = predfun,
                             col = colorPalette[1]) +
      ggplot2::stat_function(fun = function(x) predfun(x) + dtd,
                             col = colorPalette[2]) +
      ggplot2::theme_bw() +
      ggplot2::scale_colour_manual(values = colorPalette) +
      ggplot2::ylab(ylab) +
      ggplot2::xlab(xlab)
  }
  p
}


#' Plot mean comparisons for circGLM objects
#'
#' If the main predictors of interest for the circGLM are categorical, it can be
#' insightful to plot the posteriors of the group means side-by-side, which this
#' function does. This is particularly useful for ANOVA or ANCOVA type designs.
#'
#' If there are linear predictors in the model as well, the posteriors displayed
#' will correspond to the intercept parameter for each group.
#'
#' @param m A circGLM object.
#' @param alpha The alpha level of the plotted densities.
#' @param xlab The label of the x-axis.
#'
#' @export
#'
#' @seealso \code{\link{"plot_trace.circGLM"}},
#'   \code{\link{"plot_tracestack.circGLM"}},
#'   \code{\link{"plot_predict.circGLM"}},
#'   \code{\link{"plot_meanboxplot.circGLM"}},
#'   \code{\link{"plot.circGLM"}}.
#'
#' @examples
#' dat <- generateCircGLMData(nconpred = 0)
#' m   <- circGLM(th = dat[, 1], X = dat[, -1])
#' plot_meancompare.circGLM(m)
plot_meancompare.circGLM <- function(m, alpha = .7, xlab = "Mean direction") {

  ggplot2::ggplot(data = reshape2::melt(m$mu_chain, varnames = c("ID", "mu")),
                  mapping = ggplot2::aes(x = value, fill = factor(mu))) +
    ggplot2::geom_density(alpha = alpha) +
    ggplot2::scale_fill_discrete(guide = ggplot2::guide_legend(title = "Group")) +
    ggplot2::theme_bw() + ggplot2::xlab(xlab) + ggplot2::ylab("Probability")
}



#' Plot mean comparison boxplot from circGLM objects
#'
#' If the main predictors of interest for the circGLM are categorical, it can be
#' insightful to plot the posteriors of the group means side-by-side, which this
#' function does. This is particularly useful for ANOVA or ANCOVA type designs.
#'
#' If there are linear predictors in the model as well, the posteriors displayed
#' will correspond to the intercept parameter for each group.
#'
#' Some caution is needed, as a regular linear boxplot is printed, which may not
#' always be meaningful for a circular variable.
#'
#' @param m A circGLM object.
#' @param xlab The label of the x-axis.
#'
#' @export
#'
#' @seealso \code{\link{"plot_trace.circGLM"}},
#'   \code{\link{"plot_tracestack.circGLM"}},
#'   \code{\link{"plot_predict.circGLM"}},
#'   \code{\link{"plot_meancompare.circGLM"}},
#'   \code{\link{"plot.circGLM"}}.
#'
#' @examples
#' dat <- generateCircGLMData(nconpred = 0)
#' m   <- circGLM(th = dat[, 1], X = dat[, -1])
#' plot_meancompare.circGLM(m)
plot_meanboxplot.circGLM <- function(m, xlab = "Mean direction") {

  df <- reshape2::melt(m$mu_chain, varnames = c("ID", "mu"))

  ggplot2::ggplot(data = df,
                  mapping = ggplot2::aes(x     = factor(mu),
                                         y     = value,
                                         fill  = factor(mu))) +
    ggplot2::geom_boxplot(outlier.size = 0) + ggplot2::coord_flip() +
    ggplot2::scale_fill_discrete(guide = ggplot2::guide_legend(title = "Group")) +
    ggplot2::theme_bw() + ggplot2::xlab("") + ggplot2::ylab(xlab)
}


#' Make traceplots for circGLM
#'
#'  Plot traceplots from a circGLM object. This plotting method uses the standard \code{coda} traceplots.
#'
#' @param m A circGLM object.
#' @param params An optional character vector containing the parameter chains to display. If left empty, all are plotted.
#' @param ... Additional parameters passed to \code{\link[coda]{plot.mcmc}} from the coda package.
#'
#' @export
#'
#' @seealso \code{\link{"plot_tracestack.circGLM"}},
#'   \code{\link{"plot_predict.circGLM"}},
#'   \code{\link{"plot_meancompare.circGLM"}},
#'   \code{\link{"plot_meanboxplot.circGLM"}},
#'   \code{\link{"plot.circGLM"}}.
#'
#' @examples
#' plot_trace.circGLM(circGLM(rvmc(10, 1, 1)))
#'
#' dat <- generateCircGLMData()
#' plot(circGLM(th = dat[, 1], X = dat[, -1]), type = "trace")
#'
plot_trace.circGLM <- function(m, params, ...) {
  if (missing(params)) {
    plot(m$all_chains, ...)
  } else {
    plot(m$all_chains[, params], ...)
  }
  invisible(NULL)
}


#' An alternative option to plot traceplots from circGLM objects.
#'
#' @param m A circGLM object.
#' @param coef A character string, either "Beta" or "Zeta", determining whether
#'   the continuous regression predictors are shown in reparametrized form or
#'   not.
#' @param labelFormat A character vector, either "default", numbered" or
#'   "latex". By default, we find the names of the variables in the circGLM
#'   object. If \code{"numbered"}, the parameter names are numberd. The "latex"
#'   labels are useful if \code{knitr} is used with a Tikz device.
#' @param ggTheme The ggplot theme to use if
#' @param res The maximum number iterations to print. If \code{res} is larger
#'   than the number of iterations in the circGLM object, a subset of size
#'   \code{res} is selected, and it is attempted to equally space the selected
#'   iterations from the full set. This is useful if there is a very large
#'   posterior sample due to having very little thinning.
#' @param burnThinLabel Logical; if TRUE, the x-label will reflect the fact that
#'   a burn-in and a thinning factor were used. If FALSE, the x-labels will run
#'   from 1 to Q.
#'
#' @return A ggplot2 plot.
#' @export
#'
#'
#' @seealso \code{\link{"plot_trace.circGLM"}},
#'   \code{\link{"plot_predict.circGLM"}},
#'   \code{\link{"plot_meancompare.circGLM"}},
#'   \code{\link{"plot_meanboxplot.circGLM"}}, \code{\link{"plot.circGLM"}}.
#'
#'
#'
#'
#' @examples
#' plot(circGLM(rvmc(100, 0, 1)), type = "tracestack")
#'
#' dat <- generateCircGLMData()
#' plot(circGLM(th = dat[, 1], X = dat[, -1]), type = "tracestack")
#'
plot_tracestack.circGLM <- function(m,
                                    coef="Beta",
                                    labelFormat = "default",
                                    ggTheme = ggplot2::theme_bw(),
                                    res = 10000,
                                    burnThinLabel = TRUE) {

  if (is.null(m$b0_chain)) stop("No posterior sample saved for this result.")

  if (coef == "Beta") {
    pd_chain <- m$bt_chain
  } else if (coef == "Zeta") {
    pd_chain <- m$zt_chain
  } else {
    stop("Coef type not found")
  }

  ndt <- ncol(m$dt_chain)
  npd <- ncol(pd_chain)

  # Create labels
  if (labelFormat == "default") {
    yLabB0    <- "Intercept"
    yLabKp    <- "Kappa"
    yLabDt    <- if (ndt > 0) {colnames(m$dt_meandir)} else {NULL}
    yLabPd    <- if (npd > 0) {colnames(m$bt_mean)   } else {NULL}
  } else if (labelFormat == "numbered") {
    yLabB0    <- "Beta_0"
    yLabKp    <- "Kappa"
    yLabDt    <- if (ndt > 0) {paste("Delta", 1:ndt)} else {NULL}
    yLabPd    <- if (npd > 0) {paste(coef, 1:npd)   } else {NULL}
  } else if (labelFormat == "latex") {
    yLabB0    <- "$\\beta_0$"
    yLabKp    <- "$\\kappa$"
    yLabDt    <- if (ndt > 0) {paste0("$\\delta_", 1:ndt, "$")             } else {NULL}
    yLabPd    <- if (npd > 0) {paste0("\\", tolower(coef), "_", 1:npd, "$")} else {NULL}
  } else {
    stop("Unknown label format.")
  }

  allChains <- cbind(m$b0_chain, m$kp_chain, m$dt_chain, pd_chain)
  yLabs     <- c(yLabB0, yLabKp, yLabDt, yLabPd)
  colnames(allChains) <- yLabs


  # Chain length
  Q <- m$SavedIts

  # Strip values if necessary, because plotting is slow with too many values.
  if (Q > res) {
    # Take only res values to keep plotting fast.
    idx <- round(seq(1, Q, Q/res))
    allChains <- allChains[idx, ]
  } else {
    idx <- 1:Q
  }

  # Reflect the burnin and thinning in the x-label if required.
  if (burnThinLabel) {
    scaleGeom <- ggplot2::scale_x_continuous(labels = function(x) m$burnin + x*m$thin)
  } else {
    scaleGeom <- NULL
  }

  # Stack chains on top of each other.
  longChains <- reshape2::melt(allChains, varnames = c("idx", "iPar"))
  longChains[, "idx"] <- idx

  p <- ggplot2::ggplot(data = longChains,
                       mapping = ggplot2::aes(x = idx,
                                              y = value,
                                              group = iPar)) +
    ggplot2::geom_line(col = rgb(0, 0, 0, .9)) +
    ggplot2::xlab("Iteration") +
    ggplot2::ylab("") +
    ggTheme +
    scaleGeom +
    ggplot2::coord_cartesian(xlim = c(0, Q), expand = FALSE) +
    ggplot2::facet_grid(iPar ~ ., scales = "free")

  p

}
