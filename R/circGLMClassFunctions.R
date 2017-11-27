#' Extract circGLM Coefficients
#'
#' Create a table of coefficient results from a \code{circGLM} object.
#'
#' @param object A \code{circGLM} object.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A table of coefficients with their corresponding lower and upper bounds.
#' @export
#'
#' @method coef circGLM
#'
#' @examples
#' coef(circGLM(rvmc(10, 0, 1)))
#'
coef.circGLM <- coefficients.circGLM <- function(object, ...) {
  mcmcSDs <- summary(object$all_chains)$statistics[, "SD"]

  b0 <- c(object$b0_meandir, mcmcSDs["b0_chain"], object$b0_CCI)
  kp <- c(object$kp_mode, mcmcSDs["kp_chain"], object$kp_HDI)
  bt <- cbind(t(object$bt_mean), mcmcSDs[grep("bt_chain", names(mcmcSDs))], t(object$bt_CCI))
  dt <- cbind(t(object$dt_meandir), mcmcSDs[grep("dt_chain", names(mcmcSDs))], t(object$dt_CCI))

  coefmat <- rbind(Intercept = b0, Kappa = kp, bt, dt)
  colnames(coefmat) <- c("Estimate", "SD", "LB", "UB")
  coefmat
}


#' Obtain posterior model probabilities
#'
#' Compute posterior model probabilities from odds \code{x} and a prior odds.
#'
#' @param x A vector of odds for which to obtain the posterior model probabilities.
#' @param prior_odds The prior odds.
#'
#' @return A matrix with two columns, giving the relative probabilities of the first hypothesis versus the second hypothesis.
#' @export
#'
#' @examples
#' getPMP(3)
#'
getPMP <- function(x, prior_odds = 1) {
  posterior_odds <- prior_odds * x
  cbind(posterior_odds/(1 + posterior_odds), 1/(1 + posterior_odds))
}


#' Obtain Bayes Factors from circGLM objects
#'
#' Extracts the Bayes Factors from a \code{circGLM} object..
#'
#' @param m A \code{circGLM} object.
#'
#' @return A list of tables of Bayes Factors and posterior model probabilities, where applicable.
#' @export
#'
#' @examples
#' dat <- generateCircGLMData(truebeta = c(0, .2), truedelta = c(.4, .01))
#' m   <- circGLM(th = dat[, 1], X = dat[, -1])
#' BF.circGLM(m)
#'
#' dat <- generateCircGLMData(nconpred = 0)
#' m   <- circGLM(th = dat[, 1], X = dat[, -1])
#' BF.circGLM(m)
#'
#' dat <- generateCircGLMData(ncatpred = 0)
#' m   <- circGLM(th = dat[, 1], X = dat[, -1])
#' BF.circGLM(m)
#'
BF.circGLM <- function(m, prior_odds = 1) {

  # Compute posterior model probabilities
  PMP_Beta_Ineq           <- getPMP(m$BetaBayesFactors[, 1, drop = FALSE], prior_odds = prior_odds)
  PMP_Beta_Eq             <- getPMP(m$BetaBayesFactors[, 2, drop = FALSE], prior_odds = prior_odds)
  colnames(PMP_Beta_Ineq) <- c("P(bt>0)", "P(bt<0)")
  colnames(PMP_Beta_Eq)   <- c("P(bt==0)", "P(bt=/=0)")

  PMP_Mean_Ineq           <- getPMP(m$MuBayesFactors[, 1, drop = FALSE], prior_odds = prior_odds)
  PMP_Mean_Eq             <- getPMP(m$MuBayesFactors[, 2, drop = FALSE], prior_odds = prior_odds)
  colnames(PMP_Mean_Ineq) <- c("P(mu_a>mu_b)", "P(mu_a<mu_b)")
  colnames(PMP_Mean_Eq)   <- c("P(mu_a==mu_b)", "P(mu_a, mu_b)")

  list(BF_Beta       = m$BetaBayesFactors,
       PMP_Beta_Ineq = PMP_Beta_Ineq,
       PMP_Beta_Eq   = PMP_Beta_Eq,
       BF_Mean       = m$MuBayesFactors,
       PMP_Mean_Ineq = PMP_Mean_Ineq,
       PMP_Mean_Eq   = PMP_Mean_Eq)
}


#' Obtain residuals from a circGLM object
#'
#' Computes the residuals either by taking the arc distance or the cosine
#' distance between the predictions and the observed outcomes.
#'
#' @param object A \code{circGLM} object.
#' @param type Either \code{"arc"} or \code{"cosine"}, the type of distance to
#'   take.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A numeric vector of residuals. If type is \code{"arc"}, these are
#'   angles in radians. If type is \code{"cosine"}, these are numeric values
#'   between 0 and 2.
#' @export
#'
#' @method residuals circGLM
#'
#' @examples
#' m <- circGLM(rvmc(10, 0, 1))
#' residuals(m)
#'
#' # Cosine residuals
#' residuals(m, type = "cosine")
#'
residuals.circGLM <- function(object, type = "arc", ...) {

  diffs <- abs(object$data_th - object$th_hat)

  if (type == "arc") {
    return(pmin(diffs, 2*pi - diffs))
  } else if (type == "cosine") {
    return(1 - cos(diffs))
  } else {
    stop("Unknown distance type.")
  }

}


#' Obtain a prediction function from a circGLM object
#'
#' @param object A \code{circGLM} object.
#' @param linkfun A link function to use in the analysis. Should be the same as the link function.
#'
#' @return A function that takes \code{newdata} as an argument, which must be a
#'   data frame with predictors. The predictors must be the same as used in the
#'   \code{circGLM} object and must have the same column names.
#' @export
#'
#' @examples
#' dat <- generateCircGLMData()
#' m   <- circGLM(th = dat[, 1], X = dat[, -1])
#' predfun <- predict_function.circGLM(m)
#' newd <- generateCircGLMData()
#'
#' # Predicted values of the new data.
#' predfun(newd)
#'
predict_function.circGLM <- function(object, linkfun = function(x) atanLF(x, 2) ) {

  function(newdata) {

    # Check whether there are categorical and continuous predictors.
    if (length(object$dt_meandir) == 0) {
      dpart <- 0
    } else {
      d <- newdata[, colnames(object$dt_meandir)]
      dpart <- d %*% t(object$dt_meandir)
    }
    if (length(object$bt_mean) == 0) {
      xpart <- 0
    } else {
      x <- newdata[, colnames(object$bt_mean)]
      xpart <- linkfun(x %*% t(object$bt_mean))
    }
    unname(object$b0_meandir + xpart + dpart)
  }
}


#' Obtain predictions for the circGLM model
#'
#' Obtain predictions from the original dataset, or the predictions from the fitted model on a new dataset \code{newdata}.
#'
#' @param object A \code{circGLM} object.
#' @param newdata A data frame with predictors. The predictors must be the same
#'   as used in the \code{circGLM} object and must have the same column names.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return A numeric vector with predictions.
#' @export
#'
#' @method predict circGLM
#'
#' @examples
#' dat <- generateCircGLMData()
#' m   <- circGLM(dat[, 1], X = dat[, -1])
#'
#' # Predictions for the original outcome angles.
#' predict(m)
#'
#' # Predictions for new data
#' dat2  <- generateCircGLMData()
#' predict(m, newdata = dat2)
#'
predict.circGLM <- function(object, newdata, ...) {
  if (missing(newdata)) {
    return(object$th_hat)
  } else {
    predfun <- predict_function.circGLM(object)
    predfun(newdata)
  }
}


#' Compare the information criteria of several circGLM models.
#'
#' @param ... The \code{circGLM} objects to be compared.
#' @param ICs A character vector of ICs to display.
#'
#' @return A matrix with a column of information criteria for each model.
#' @export
#'
#' @examples
#' Xcat <- c(rep(0, 5), rep(1, 5))
#' th <- rvmc(10, 0, 4) + Xcat
#'
#' # Compare a model that includes group differences with a model that does not.
#' IC_compare.circGLM(circGLM(th), circGLM(th, X = Xcat))
IC_compare.circGLM <- function(...,
                               ICs = c("n_par", "lppd",
                                       "AIC_Bayes", "DIC", "DIC_alt",
                                       "WAIC1", "WAIC2",
                                       "p_DIC", "p_DIC_alt",
                                       "p_WAIC1", "p_WAIC2")) {
  ms <- list(...)

  comtab <- sapply(ms, function(m) m[ICs])

  colnames(comtab) <- as.character(match.call())[2:(ncol(comtab) + 1)]
  comtab
}


#' Compute the median direction
#'
#' This function computes the median direction, which is defined as the middle observation of the shortest arc containing all observations.
#'
#' @param th A vector of angles in radians.
#' @param fastMethod Logical; If \code{TRUE}, the data is rotated so that the
#'   mean is \code{pi} and linear methods are applied. If \code{FALSE}, the arcs
#'   between each set of data points must be computed, which is much slower. For
#'   data that is very strongly spread out, the fast method might not give the
#'   correct value.
#'
#' @return An angle in radians, the median direction.
#' @export
#'
#' @examples
#' medianDirection(rvmc(30, 0, 2))
#'
medianDirection <- function(th, fastMethod = TRUE) {

  if (fastMethod) {

    # The fast method uses the C++ method circQuantile.
    return(as.numeric(circQuantile(th, .5)))
  } else {

    # The slower method computes the mean arc distance between all angles th.
    meandiffs <- sapply(th, function(x) {
      pi - mean(abs(pi - abs(x - th)))
    })
    return(th[which.min(meandiffs)])
  }
}



#' Estimate the modal direction
#'
#' Estimates the mode as the midpoint of the highest density interval.
#'
#' The highest density interval is computed as the shortest interval containing
#' \code{modebw} of the values in \code{th}. For circular data however, this
#' definition is not useful, and we should instead look for the shortest arc
#' that contains \code{modebw} of the data. This is done by rotating the data
#' such that the mean direction is \code{pi}, and then applying the usual linear
#' methods.
#'
#' @param th A vector of angles in radians.
#' @param modebw Numeric between 0 and 1. The modes are estimated by taking the
#'   midpoint of a highest density interval. Specifically, the mode is the
#'   midpoint of the interval that contains \code{modebw} of the values of
#'   \code{th}. Reasonable values are roughly between .005 and .2, although
#'   lower values may be reasonable there are a lot of observations in
#'   \code{th}.
#'
#' @return An angle in radians.
#' @export
#'
#' @examples
#' modalDirection(rvmc(30, 0, 2))
#'
modalDirection <- function(th, modebw = .1) {

  mdir <- computeMeanDirection(th)

  mdir - pi + estimateMode(th - mdir + pi , cip = modebw)
}



#' Compute the Circular Standard Deviation
#'
#' Returns the circular standard deviation of a vector of circular data which is
#' defined as the square root of minus 2 times the log of the mean resultant
#' length.
#'
#' @param x A vector of angles.
#'
#' @return A numeric, the circular standard deviation.
#' @export
#'
circSD <- function(x) {
  sqrt(-2 * log(sqrt(sum(cos(x))^2 + sum(sin(x))^2) / length(x)))
}


#' Obtain different central tendencies and CIs from a circGLM object
#'
#' Computes the mean (arithmetic or mean direction), median, and mode estimate
#' for the MCMC chains of a \code{circGLM} object, as well as a credible interval.
#'
#' The summary statistics computed have to be computed differently for linear
#' and circular variables.
#'
#' @param m A \code{circGLM} object.
#' @param modebw Numeric between 0 and 1. The modes are estimated by taking the
#'   midpoint of a highest density interval. Specifically, the mode is the
#'   midpoint of the interval that contains \code{modebw} of the density of the
#'   posterior. Reasonable values are roughly between .005 and .2, although
#'   lower values may be reasonable if the number of iterations, Q, is large.
#' @param ciperc The confidence interval percentage.
#'
#' @return A matrix with the parameters as rows, and on the columns central
#'   tendencies and appropriate credible intervals (circular quantiles and
#'   Highest Density Intervals).
#' @export
#'
#' @examples
#' dat <- generateCircGLMData()
#' m   <- circGLM(dat[, 1], X = dat[, -1])
#' mcmc_summary.circGLM(m)
#'
mcmc_summary.circGLM <- function(m, modebw = .1, ciperc = .95) {

  nms <- colnames(m$all_chains)

  # Obtain circular central tendencies.
  circVars   <- grep("dt|b0|mu", nms)
  circChains <- m$all_chains[, circVars, drop = FALSE]
  circCTs    <- t(apply(circChains, 2, function(x) {
    c(Mean   = computeMeanDirection(x),
      Median = medianDirection(x),
      Mode   = estimateModeCirc(x, modebw),
      SD     = circSD(x),
      computeHDICirc(x, ciperc))
  }))

  lineVars <- grep("kp|bt", nms)
  lineChains <- m$all_chains[, lineVars, drop = FALSE]
  lineCTs    <- t(apply(lineChains, 2, function(x) {
    c(Mean   = mean(x),
      Median = median(x),
      Mode   = estimateMode(x, modebw),
      SD     = sd(x),
      computeHDI(x, cip = ciperc))
  }))

  out <- rbind(circCTs, lineCTs)
  colnames(out)[5:6] <- c("LB", "UB")
  out
}



