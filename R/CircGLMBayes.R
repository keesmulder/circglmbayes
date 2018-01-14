#' circglmbayes: A package for the Bayesian circular GLM.
#'
#' This package contains functions to perform a Bayesian circular GLM, which
#' allows regressing a circular outcome on linear and categorical predictors.
#' The model used in this package is similar to the model used by
#' \code{lm.circular} form the package \code{circular}. Differences are that the model used
#' by this package treats categorical variables specially. In addition, several
#' hypothesis testing options are provided.
#'
#' Estimation and uncertainty intervals are all performed in a Bayesian manner
#' through MCMC. Bayesian hypothesis tests are provided through the Bayes
#' factor.
#'
#' @section Functions: The main function of the package is
#'   \code{\link{circGLM}}, which runs an MCMC sampler in \code{C++} through
#'   \code{Rcpp}. This sampler returns an S3 object of type
#'   \code{circGLM}, which can be further analyzed through associated
#'   \code{\link{plot.circGLM}} and \code{\link{print.circGLM}} functions.
#'
#' @useDynLib circglmbayes, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom grDevices rgb
#' @importFrom graphics plot
#' @importFrom stats coef median rchisq rnorm runif sd model.matrix
#' @docType package
#' @name circglmbayes
NULL



.onUnload <- function (libpath) {
  library.dynam.unload("circglmbayes", libpath)
}
