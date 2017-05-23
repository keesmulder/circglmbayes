#' Check if a predictor is dichotomous
#'
#' @param x A character or numerical vector to be tested.
#'
#' @return A logical, TRUE if the \code{x} has dummy coding (0, 1), FALSE otherwise.
#'
is.dichotomous <- function(x) {
  if (length(unique(x)) == 2) {
    if (all(x == 0 | x == 1)) {
      return(TRUE)
    } else {
      warning("A predictor might be dichotomous but not 0|1.")
    }
  }
  FALSE
}


#' Fix names for circGLM vector output
#'
#' A function to change the names produced in the Rcpp code to more human readable forms.
#'
#' This is only done if the \code{\link{circGLM}} function is used with \code{output = "vector"}.
#'
#' @param nms The original names.
#'
#' @return A character vector of the same length as \code{nms}.
fixResultNames <- function(nms){

  nms[grep("b0_CCI", nms)]     <- c("b0_CCI_LB", "b0_CCI_UB")
  nms[grep("kp_HDI", nms)]     <- c("kp_HDI_LB", "kp_HDI_UB")

  # Beta/Zeta
  nbts <- length(grep("bt_mean", nms))
  if (nbts > 0) {
    nms[grep("bt_mean", nms)]    <- paste0("bt_", 1:nbts, "_mean")
    nms[grep("zt_mean", nms)]    <- paste0("zt_", 1:nbts, "_mean")
    nms[grep("zt_mdir", nms)]    <- paste0("zt_", 1:nbts, "_mdir")
    nms[grep("bt_propacc", nms)] <- paste0("bt_", 1:nbts, "_propacc")
    nms[grep("bt_CCI", nms)] <- paste0("bt_",
                                       rep(1:nbts, each=2),
                                       c("_LB", "_UB"))
    nms[grep("zt_CCI", nms)] <- paste0("zt_",
                                       rep(1:nbts, each=2),
                                       c("_LB", "_UB"))
  }

  # Delta
  ndts <- length(grep("dt_meandir", nms))
  if (ndts > 0) {

    nms[grep("dt_meandir", nms)] <- paste0("dt_", 1:ndts, "_mdir")
    nms[grep("dt_propacc", nms)] <- paste0("dt_", 1:ndts, "_propacc")
    nms[grep("dt_CCI", nms)] <- paste0("dt_",
                                       rep(1:ndts, each = 2),
                                       c("_LB", "_UB"))
  }

  nms
}


#' Fitting Bayesian circular outcome General Linear Models
#'
#' The main function for running Bayesian circular GLMs.
#'
#' \code{circGLM} performs an mcmc sampler that generates a sample from the
#' posterior of the intercept \eqn{\beta_0}, regression coefficients
#' \eqn{\beta}, group mean direction differences \eqn{\delta} and residual
#' \eqn{\kappa}.
#'
#' An attempt is made to split the predictor matrix \code{X} into continuous and
#' categorical predictors. This is done so that the categorical predictors can
#' be treated differently, which removes the arbitrary dependence on the
#' labeling of categorical predictors and ensures that each group has a
#' regression line of the same shape.
#'
#' The main results obtained are estimates and credible intervals for the
#' parameters, posterior samples, and Bayes factors for various standard
#' hypothesis comparisons.
#'
#' As with all mcmc samplers, convergence must be checked, and tuning parameters
#' \code{bwb} and \code{reparametrize} can be tweaked if the sampler converges
#' badly. The circGLM object that is returned contains proportions acceptepted
#' which can be used to monitor performance.
#'
#'
#'
#' @param th A vector of angles in radians or degrees, representing the circular
#'   outcome we want to predict. If any value is larger than \code{2 * pi}, the
#'   input is transformed to radians. Otherwise, \code{th} is treated as
#'   radians.
#' @param X A matrix of predictors, both continuous (linear) and categorical.
#'   Categorical predictors must be in (0, 1), so that they are treated
#'   separately from the continuous predictors. If not, or if
#'   \cite{skipDichSplit = TRUE}, they will be treated as linear predictors.
#'   continuous.
#' @param conj_prior A numeric vector of length 3, containing, in that order,
#'   prior mean direction, prior resultant length, and prior sample size. Used
#'   for the von Mises part of the model, beta_0 and kappa.
#' @param bt_prior_musd A numeric vector of length 2, or \code{NA}. If
#'   \code{bt_prior_musd = NA}, a constant prior is used. If it is a numeric
#'   vector of length 2, the first value is the mean, and the second value is
#'   the standard deviation.
#' @param starting_values A numeric vector with starting values for the mcmc
#'   sampler. The length of the numeric vector should be 2 plus the number of
#'   columns in X.
#' @param bwb A numeric vector, where the length is at least the number of
#'   continuous predictors. This is a tuning parameters used in sampling of
#'   beta. New values are sampled uniformly around the current value of beta
#'   with bounds at \code{bt_cur - bwb} and \code{bt_cur + bwb}. If
#'   \code{reparametrize = TRUE}, bwb corresponds to the bounds around the
#'   reparametrized values.
#' @param Q Integer; The number of iterations to perform.
#' @param burnin Integer; The number of burnin (warmup) iterations.
#' @param thin Integer; The number of parameters sets to sample for each
#'   parameter set that is saved. Can be used to save memory if \code{Q} is
#'   large.
#' @param kappaModeEstBandwith Numeric between 0 and 1. The mode of \code{kappa}
#'   is estimated by taking the midpoint of a highest density interval.
#'   Specifically, it is the midpoint of the interval that contains
#'   \code{kappaModeEstBandwith} of the density of the posterior. Reasonable
#'   values are roughly between .005 and .2, although lower values may be
#'   reasonable if Q is large.
#' @param CIsize The size of the credible intervals. This is used for all
#'   parameters, whether they use highest density intervals, circular quantiles
#'   or regular quantiles.
#' @param r A numeric. \code{r} is the parameter used in the link function
#'   \eqn{g(x, r) = r atan(x)}. If \code{r = 2}, the link function maps the real
#'   line to the full circle. If \code{r < 2} the link functions maps to a
#'   proportion \code{r / 2} of the circle. If \code{r > 2}, the link functions
#'   can reach the same are of the circle multiple times, which is unlikely to
#'   be useful, and should be used with caution.
#' @param returnPostSample Logical indicating whether the mcmc sample itself
#'   should be returned. Should only be set to FALSE if there are memory
#'   constraints, as many subsequent analyses rely on the posterior sample
#'   directly.
#' @param output A character string, either \code{"list"} or \code{"vector"}. In
#'   most situations, \code{"list"} should be used, which returns a circGLM
#'   object. The \code{"vector"} options is only useful for simulation studies
#'   etc.
#' @param reparametrize Logical; If \code{TRUE}, proposals for beta are drawn
#'   uniformly around a reparametrization \code{zt = pi * atan(bt) / 2}, so from
#'   \code{zt_can = runif(1, zt - bwb, zt + bwb)}, which is then transformed
#'   back. Then, the proposals amount to the truncated cauchy pdf. If
#'   \code{FALSE}, proposals for beta are drawn on uniformly around beta, so
#'   from \code{bt_can = runif(1, bt_cur - bwb, bt_cur + bwb)}.
#' @param groupMeanComparisons Logical indicating whether mean comparisons in
#'   the form of Bayes Factors and posterior model probabilites should be
#'   computed.
#' @param skipDichSplit Logical indicating whether to treat categorical
#'   predictor specially. Usually, \code{skipDichSplit = TRUE} should be used.
#'   In practice, this removes the arbitrary dependence on the labeling of
#'   categorical predictors and ensures that each group has a regression line of
#'   the same shape.
#' @param centerOnly Logical; If \code{TRUE}, the continuous predictors are
#'   centered only, not standardized. If \code{FALSE}, the continuous predictors
#'   are standardized.
#'
#' @return A circGLM object, which can be further analyzed with its associated
#'   \code{\link[plot.circGLM]{plot}} and \code{\link[print.circGLM]{print}}
#'   functions.
#' @export
#'
#' @seealso \code{\link{"print.circGLM"}}, \code{\link{"plot.circGLM"}},
#'   \code{\link{"coef.circGLM"}}, \code{\link{"BF.circGLM"}},
#'   \code{\link{"residuals.circGLM"}}, \code{\link{"predict.circGLM"}},
#'   \code{\link{"predict_function.circGLM"}},
#'   \code{\link{"mcmc_summary.circGLM"}}, \code{\link{"IC_compare.circGLM"}}.
#'
#' @examples
#' dat <- generateCircGLMData()
#' m   <- circGLM(th = dat[, 1], X = dat[, -1])
#' print(m)
#' print(m, type = "all")
#' plot(m, type = "tracestack")
circGLM <- function(th,
                    X = matrix(nrow = length(th), ncol = 0),
                    conj_prior = rep(0, 3),
                    bt_prior_musd = c("mu" = 0, "sd" = 1),
                    starting_values = c(0, 1, rep(0, ncol(X))),
                    bwb = rep(.05, ncol(X)),
                    Q = 10000,
                    burnin = 1000,
                    thin = 1,
                    kappaModeEstBandwith = .1,
                    CIsize = .95,
                    r = 2,
                    returnPostSample = TRUE,
                    output = "list",
                    reparametrize = FALSE,
                    groupMeanComparisons = TRUE,
                    skipDichSplit = FALSE,
                    centerOnly = FALSE) {

  # Check if the inputs are matrices.
  if (!is.matrix(th)) th <- as.matrix(th)
  if (!is.matrix(X))  X <- as.matrix(X)

  # Check if theta is in radians
  if (any(th > 2*pi)) {
    th <- th * pi / 180
  }

  # Indices of dichotomous predictors.
  dichInd <- apply(X, 2, is.dichotomous)

  if (!skipDichSplit) {

    # Standardize continuous predictors.
    stX <- scale(X[, !dichInd, drop = FALSE], scale = !centerOnly)

    d <- X[, dichInd, drop = FALSE]

  } else {
    # Standardize continuous predictors.
    stX <- cbind(scale(X[, !dichInd, drop = FALSE], scale = !centerOnly),
                 X[, dichInd, drop = FALSE])

    d <- matrix(nrow = nrow(X), ncol = 0)
  }

  # Names of the continuous and categorical predictors.
  nbt <- ncol(stX)
  ndt <- ncol(d)
  if (is.null(colnames(stX))) {
    bt_names <- paste0("bt_", 1:nbt)
    zt_names <- paste0("zt_", 1:nbt)
  } else {
    bt_names <- colnames(stX)
    zt_names <- paste0(colnames(stX), "_zt")
  }
  if (is.null(colnames(d))) {
    dt_names <- paste0("dt_", 1:ndt)
    mu_names <- c("Reference", paste0("mu_", 1:ndt))
  } else {
    dt_names <- colnames(d)
    mu_names <- c("Reference", dt_names)
  }

  # Set the a prior for beta if needed.
  if (ncol(stX) > 0) {
    bt_prior <- matrix(bt_prior_musd, nrow = ncol(stX), ncol = 2, byrow = TRUE)
  } else {
    bt_prior <- matrix(NA, nrow = 0, ncol = 2, byrow = TRUE)
  }

  # Run the C++ mcmc sampler.
  res <- circGLMC(th = th, X = stX, D = d,
                  conj_prior = conj_prior,
                  bt_prior = bt_prior,
                  starting_values = starting_values,
                  burnin = burnin, thin = thin, bwb = bwb,
                  kappaModeEstBandwith = kappaModeEstBandwith,
                  CIsize = CIsize,
                  Q = Q, r = r,
                  returnPostSample = returnPostSample,
                  bt_prior_type = as.numeric(!is.na(bt_prior_musd)[1]),
                  reparametrize = reparametrize,
                  groupMeanComparisons = groupMeanComparisons)


  ### FIXING NAMES

  # Set some names for clarity in the output.
  colnames(res$b0_CCI)     <- "Beta_0"
  rownames(res$b0_CCI)     <- c("LB", "UB")
  colnames(res$kp_HDI)     <- "Kappa"
  rownames(res$kp_HDI)     <- c("LB", "UB")

  # Set names for beta only if there are beta's.
  if (length(res$bt_mean) > 0) {
    colnames(res$bt_propacc) <- colnames(res$bt_CCI) <-
      colnames(res$bt_mean) <- bt_names
    rownames(res$bt_CCI) <- rownames(res$zt_CCI) <- c("LB", "UB")

    if (returnPostSample) colnames(res$bt_chain) <- bt_names

    colnames(res$zt_CCI) <- colnames(res$zt_mdir) <- colnames(res$zt_mean) <- zt_names

    # Fix names for Beta Bayes Factors
    rownames(res$BetaIneqBayesFactors) <- rownames(res$BetaSDDBayesFactors) <- bt_names
    res$BetaBayesFactors <- cbind(res$BetaIneqBayesFactors, res$BetaSDDBayesFactors)
    colnames(res$BetaBayesFactors) <- c("BF(bt>0:bt<0)",
                                        "BF(bt==0:bt=/=0)")

  }



  # Set names for delta only if there are delta's.
  if (length(res$dt_meandir) > 0) {

    # Fix names for delta estimates
    colnames(res$dt_meandir) <- rownames(res$DeltaIneqBayesFactors) <-
      colnames(res$dt_propacc) <- colnames(res$dt_CCI) <- dt_names

    if (returnPostSample) colnames(res$dt_chain) <- dt_names

    rownames(res$dt_meandir) <- "MeanDir"
    rownames(res$dt_CCI)  <- c("LB", "UB")
    rownames(res$dt_propacc) <- "ProportionAccepted"

    if (returnPostSample) colnames(res$mu_chain) <- mu_names

    # Fix names for Delta Ineq Bayes Factors
    colnames(res$DeltaIneqBayesFactors) <- "BF(dt>0:dt<0)"

    if (groupMeanComparisons) {

      ngroup  <- sum(dichInd) + 1
      basemat <- matrix(1:ngroup, ncol = ngroup, nrow = ngroup)
      first   <- t(basemat)[lower.tri(basemat)]
      last    <- basemat[lower.tri(basemat, diag = FALSE)]

      rownames(res$MuIneqBayesFactors) <- rownames(res$MuSDDBayesFactors) <-
        paste0("[", mu_names[first], ", ", mu_names[last],"]")
      res$MuBayesFactors <- cbind(res$MuIneqBayesFactors,
                                  res$MuSDDBayesFactors)
      colnames(res$MuBayesFactors) <- c("BF(mu_a>mu_b:mu_a<mu_b)", "BF(mu_a==mu_b:(mu_a, mu_b))")

      names(dimnames(res$MuBayesFactors)) <- c("[mu_a, mu_b]", "Comparison")
    }
  }

  rownames(res$TimeTaken) <- c("Initialization", "Loop", "Post-processing", "Total")
  colnames(res$TimeTaken) <- "Time (sec)"

  # Define a coda mcmc object containing all chains.
  res$all_chains <- coda::mcmc(as.data.frame(res[grep("chain", names(res))]),
                         start = burnin, thin = thin)

  # Choose how to return the output.
  if (output == "list") {

    # Add a class 'circGLM', so that we can use print and plot methods for this class.
    class(res) <- c("circGLM", class(res))


    # Save some of the inputs in the output for easy access.
    res$Call      <- match.call()
    res$thin      <- thin
    res$burnin    <- burnin
    res$data_th   <- th
    res$data_X    <- X
    res$data_d    <- d
    res$data_stX  <- stX
    res$r         <- r

    return(res)

  } else if (output == "vector") {
    if (returnPostSample == TRUE) {message("Vector output with full chains.")}

    out <- unlist(res)
    names(out) <- fixResultNames(names(out))
    return(out)

  } else {
    stop(paste("Output type", output, "not found"))
  }
}
