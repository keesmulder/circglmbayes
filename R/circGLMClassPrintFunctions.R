#' Print circGLM Object
#'
#' General print function for \code{circGLM} objects, which dispatches the chosen type
#' of printing to the corresponding function.
#'
#' @param x A \code{circGLM} object to be printed.
#' @param type Character string giving the type of printing, such as
#'   \code{"text"}, \code{"mcmc"}, \code{"all"}, \code{"coef"}.
#' @param ... Additional arguments to be passed to print functions.
#'
#' @export
#'
#' @seealso \code{\link{print_text.circGLM}},
#'   \code{\link{print_mcmc.circGLM}},
#'   \code{\link{print_all.circGLM}},
#'   \code{\link{print_coef.circGLM}}.
#'
#'
#' @examples
#' print(circGLM(th = rvmc(10, 1, 1)))
#'
#' dat <- generateCircGLMData()
#' cglmmod <- circGLM(th ~ ., dat)
#'
#' print(cglmmod)
#'
#' print(cglmmod, type = "mcmc")
#'
#' print(cglmmod, type = "all")
#'
#' print(cglmmod, type = "coef")
#'
print.circGLM <- function(x, type = "text", ...) {
  printFunName <- paste0("print_", type, ".circGLM")
  do.call(printFunName, args = c(list(m = x), list(...)))
}


#' Print the main results from a \code{circGLM} object.
#'
#' @param m A \code{circGLM} object.
#' @param digits Number of digits to display.
#'
#' @export
#'
#' @seealso \code{\link{print_mcmc.circGLM}},
#'   \code{\link{print_all.circGLM}},
#'   \code{\link{print_coef.circGLM}},
#'   \code{\link{print.circGLM}}.
#'
#' @examples
#' print(circGLM(th = rvmc(10, 1, 1)), type = "text")
#'
#' dat <- generateCircGLMData()
#' cglmmod <- circGLM(th = dat[, 1], X = dat[, -1])
#' print(cglmmod, type = "text")
print_text.circGLM <- function(m, digits = 3) {
  cat("Bayesian circular GLM \n")
  cat("\nCall:\n",
      paste(deparse(m$Call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("MCMC run for", m$TotalIts, "its,", m$SavedIts, "used. \n\n")
  cat("Coefficients:\n")
  print(round(coef(m), digits))
  cat("\nDIC: ", round(m$DIC, digits), "\n")
  cat("WAIC:", round(m$WAIC1, digits), "\n")

}


#' Print the mcmc results from a circGLM object
#'
#' This prints a number of diagnostics about the results of a \code{circGLM} objects
#' through  \code{\link[coda]{summary.mcmc}} from the \code{coda} package. In
#' particular, the standard errors may be of interest.
#'
#' Note that the standard error and convergence diagnostics computed by
#' \code{coda} are not necessarily trustworthy.
#'
#' @param m A \code{circGLM} object.
#' @param ... Additional arguments to be passed to coda printing functions.
#'
#' @export
#'
#' @seealso \code{\link{print_text.circGLM}},
#'   \code{\link{print_all.circGLM}},
#'   \code{\link{print_coef.circGLM}},
#'   \code{\link{print.circGLM}}.
#'
#' @examples
#' print(circGLM(th = rvmc(10, 1, 1)), type = "mcmc", digits = 3)
#'
#' dat <- generateCircGLMData()
#' cglmmod <- circGLM(th = dat[, 1], X = dat[, -1])
#' print(cglmmod, type = "mcmc")
print_mcmc.circGLM <- function(m, ...) {
  print(summary(m$all_chains), ...)
}


#' Print all results from a circGLM object
#'
#' This function prints the full list of results from a \code{circGLM} object. The function extracts all
#' the scalar results and displays these together, then prints all further list elements. The full
#' chains are not printed.
#'
#' @param m A \code{circGLM} object.
#' @param digits Number of digits to display.
#'
#' @export
#'
#' @seealso \code{\link{print_text.circGLM}}, \code{\link{print_mcmc.circGLM}},
#'   \code{\link{print_coef.circGLM}}, \code{\link{print.circGLM}}.
#'
#' @examples
#' print(circGLM(th = rvmc(10, 1, 1)), type = "all")
#'
#' dat <- generateCircGLMData()
#' cglmmod <- circGLM(th ~ ., dat)
#' print(cglmmod, type = "all")
print_all.circGLM <- function(m, digits = 3) {

  # Remove empty results (ie. parameters not in current model)
  empties <- lapply(m, length) == 0
  if (any(empties)) m <- m[-which(empties)]

  # Gather single results
  singles <- (lapply(m, length) == 1) & (!sapply(m, is.character))
  singmat <- as.matrix(round(unlist(m[singles]), digits))
  print(singmat)

  cat("\n\n")

  # Print the rest
  if (any(singles)) m <- m[-which(singles)]

  rem_str <- "chain|Call|data_|curpars|th_hat|MuSDDBayesFactors"
  print(m[-grep(rem_str, names(m))])

  cat("\n\n")
  return(invisible(NULL))
}


#' Print circGLM coefficients
#'
#' @param m A \code{circGLM} object.
#' @param digits Number of digits to display.
#'
#' @export
#'
#' @seealso \code{\link{print_text.circGLM}},
#'   \code{\link{print_mcmc.circGLM}},
#'   \code{\link{print_all.circGLM}},
#'   \code{\link{print.circGLM}}.
#'
#' @examples
#' print(circGLM(th = rvmc(10, 0, 1)), type = "coef")
#'
#' dat <- generateCircGLMData()
#' cglmmod <- circGLM(th = dat[, 1], X = dat[, -1])
#' print(cglmmod, type = "coef")
print_coef.circGLM <- function(m, digits = 3) {
  print(round(coef(m), digits))
}


