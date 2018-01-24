#' Generate data that follows the circular GLM model
#'
#' This function samples data according to the circular GLM model. A set of true
#' values for the parameters can be entered, and a dataset is returned that is
#' drawn from the corresponding model. The link function can also be selected.
#'
#' This function can also be used as a wrapper for sampling von Mises data, if
#' \code{nconpred = 0}, \code{ncatpred = 0}. Then, \code{beta_0} is the mean of
#' the von Mises distribution and \code{residkappa} is the concentration
#' parameter \eqn{\kappa}.
#'
#' In order to make this function more useful in simulations, the true
#' parameters are also added to the data set that is returned as attributes.
#'
#' @param n Integer; the sample size to be generated.
#' @param residkappa A non-negative numeric; the residual concentration
#'   parameter. This is the \eqn{\kappa} of the von Mises distribution that the
#'   residuals follow.
#' @param nconpred Integer; The number of continuous (linear) predictors to be
#'   generated.
#' @param ncatpred Integer; The number of categorical predictors to be
#'   generated.
#' @param truebeta0 An angle in radians representing \code{beta_0}, which
#'   functions as the intercept.
#' @param truebeta A numeric vector containing the values for the regression
#'   coefficients of the continuous predictors.
#' @param truedelta A numeric vector containing angles in radians that represent
#'   the group differences for each of the categorical predictors.
#' @param linkfun The link function to use. The default is the canonical
#'   arctangent link.
#'
#' @return A numeric matrix containing a dataset sampled according to the
#'   circular GLM model. The first column \code{th} represents the circular
#'   outcome in radians. The following columns represent the linear predictors
#'   and are named \code{l1}, \code{l2}, ... . The following columns represent
#'   the categorical predictors and are named \code{c1}, \code{c2}, ... . The
#'   matrix also has attributes containing the true values of the parameters,
#'   the used link function, and a proportion \code{u} showing the proportion of
#'   the data that is on the semicircle closest to \code{beta_0}.
#' @export
#'
#' @examples
#'
#' # Von Mises data with mean 2, kappa 3.
#' generateCircGLMData(truebeta0 = 2, residkappa = 3,
#'                     nconpred = 0, ncatpred = 0)
#'
#' # circGLM data
#' generateCircGLMData(n = 20, nconpred = 4, truebeta = c(0, 0.4, 0.2, 0.05))
#' 
generateCircGLMData <- function(n = 30, residkappa = 5,
                                nconpred = 2, ncatpred = 2,
                                truebeta0 = pi/2,
                                truebeta  = rep(.25, nconpred),
                                truedelta = rep(1, ncatpred),
                                linkfun   = function(x) 2 * atan(x)) {

  dtpart <- btpart <- 0

  Xcon <- matrix(nrow = n, ncol = nconpred)
  Xcat <- matrix(nrow = n, ncol = ncatpred)

  # Check whether true parameter values must be drawn.
  if (!is.numeric(truebeta0)  && truebeta0 == "random") {
    truebeta0 <- runif(1, -2, 8)
  }
  if (!is.numeric(truebeta)   && truebeta == "random") {
    truebeta  <- rnorm(nconpred + ncatpred, sd = 0.5)
  }
  if (!is.numeric(residkappa) && residkappa == "random") {
    residkappa  <- rchisq(1, 10)
  }

  # Generate predictors
  if (nconpred > 0) {
    Xcon <- sapply(1:nconpred, function(x) scale(rnorm(n)))
    colnames(Xcon)  <- paste0("l", 1:nconpred)
    btpart <- linkfun(apply(Xcon, 1, "%*%", truebeta))
  }

  if (ncatpred > 0) {
    Xcat  <- sapply(1:ncatpred, function(x) sample(0:1, size = n, 
                                                   replace = TRUE))
    colnames(Xcat) <- paste0("c", 1:ncatpred)
    dtpart <- apply(Xcat, 1, "%*%", truedelta)
  }



  # Generate values for the circular outcome.
  thpred <- truebeta0 + dtpart + btpart

  therr  <- rvmc(n, 0, residkappa)
  th     <- (thpred + therr) %% (2*pi)

  dmat   <- cbind(th, Xcon, Xcat)

  # Save the percentage of data that is found on the half circle closest to the
  # true beta_0.
  uLB <- truebeta0 - pi/2 %% (2*pi)
  uUB <- truebeta0 + pi/2 %% (2*pi)
  if (uLB < uUB) {
    u <- mean(uLB < th & th < uUB)
  } else {
    u <- mean(uLB < th | th < uUB)
  }

  # Add the true values as attributes.
  attr(dmat, "truebeta0")  <- truebeta0
  attr(dmat, "truebeta")   <- truebeta
  attr(dmat, "truezeta")   <- (2/pi)*atan(truebeta)
  attr(dmat, "residkappa") <- residkappa
  attr(dmat, "linkfun")    <- linkfun
  attr(dmat, "u")          <- u

  return(dmat)
}

