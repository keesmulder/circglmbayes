#' Extract circGLM Coefficients
#'
#' Create a table of coefficient results from a circGLM object.
#'
#' @param m A circGLM object.
#'
#' @return A table of coefficients with their corresponding lower and upper bounds.
#' @export
#'
#' @examples
#' coef(circGLM(rvmc(10, 0, 1)))
#'
coef.circGLM <- coefficients.circGLM <- function(m) {
  mcmcSDs <- summary(m$all_chains)$statistics[, "SD"]

  b0 <- c(m$b0_meandir, mcmcSDs["b0_chain"], m$b0_CCI)
  kp <- c(m$kp_mode, mcmcSDs["kp_chain"], m$kp_HDI)
  bt <- cbind(t(m$bt_mean), mcmcSDs[grep("bt_chain", names(mcmcSDs))], t(m$bt_CCI))
  dt <- cbind(t(m$dt_meandir), mcmcSDs[grep("dt_chain", names(mcmcSDs))], t(m$dt_CCI))

  coefmat <- rbind(Intercept = b0, Kappa = kp, bt, dt)
  colnames(coefmat) <- c("Estimate", "SD", "LB", "UB")
  coefmat
}


#' Obtain Bayes Factors from circGLM objects
#'
#' @param m A circGLM object
#' @param digits
#'
#' @return A list of tables of Bayes Factors, where applicable.
#' @export
#'
#' @examples
#' dat <- generateCircGLMData(truebeta = c(0, .2), truedelta = c(.4, .01))
#' m   <- circGLM(dat[, 1], X = dat[, -1])
#' BF.circGLM(m)
#'
#' dat <- generateCircGLMData(nconpred = 0)
#' m   <- circGLM(dat[, 1], X = dat[, -1])
#' BF.circGLM(m)
#'
#' dat <- generateCircGLMData(ncatpred = 0)
#' m   <- circGLM(dat[, 1], X = dat[, -1])
#' BF.circGLM(m)
#'
BF.circGLM <- function(m) {
  list(Beta = m$BetaBayesFactors,
       Mean = m$MuBayesFactors)
}


#' Obtain residuals from a circGLM object
#'
#' Computes the residuals either by taking the arc distance or the cosine
#' distance between the predictions and the observed outcomes.
#'
#' @param m A circGLM object.
#' @param type Either \code{"arc"} or \code{"cosine"}, the type of distance to
#'   take.
#'
#' @return A numeric vector of residuals. If type is \code{"arc"}, these are
#'   angles in radians. If type is \code{"cosine"}, these are numeric values
#'   between 0 and 2.
#' @export
#'
#' @examples
#' m <- circGLM(rvmc(10, 0, 1))
#' residuals.circGLM(m)
#' residuals.circGLM(m, type = "cosine")
#'
residuals.circGLM <- function(m, type = "arc") {

  diffs <- abs(m$data_th - m$th_hat)

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
#' @param m A circGLM object.
#'
#' @return A function that takes \code{newdata} as an argument, which must be a
#'   data frame with predictors. The predictors must be the same as used in the
#'   circGLM object and must have the same column names.
#' @export
#'
#' @examples
#' dat <- generateCircGLMData()
#' m   <- circGLM(th = dat[, 1], X = dat[, -1])
#' predict_function.circGLM(m)
predict_function.circGLM <- function(m, linkfun = function(x) atanLF(x, 2) ) {

  function(newdata) {

    # Check whether there are categorical and continuous predictors.
    if (length(m$dt_meandir) == 0) {
      dpart <- 0
    } else {
      d <- newdata[, colnames(m$dt_meandir)]
      dpart <- d %*% t(m$dt_meandir)
    }
    if (length(m$bt_mean) == 0) {
      xpart <- 0
    } else {
      x <- newdata[, colnames(m$bt_mean)]
      xpart <- linkfun(x %*% t(m$bt_mean))
    }
    m$b0_meandir + xpart + dpart
  }
}


#' Obtain predictions for the circGLM model
#'
#' @param m A circGLM object.
#' @param newdata A data frame with predictors. The predictors must be the same
#'   as used in the circGLM object and must have the same column names.
#'
#' @return A numeric vector with predictions.
#' @export
#'
#' @examples
#' dat <- generateCircGLMData()
#' m   <- circGLM(dat[, 1], X = dat[, -1])
#'
#' # Predictions for the original outcome angles.
#' predict.circGLM(m)
#'
#' # Predictions for new data
#' dat2  <- generateCircGLMData()
#' predict.circGLM(m, newdata = dat2)
predict.circGLM <- function(m, newdata) {
  if (missing(newdata)) {
    return(m$th_hat)
  } else {
    predfun <- predict_function.circGLM(m)
    predfun(newdata)
  }
}


#' Compare the information criteria of several circGLM models.
#'
#' @param ... The circGLM objects to be compared.
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

