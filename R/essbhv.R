#' Basic Human Values data
#'
#' Data from Dutch respondents of the European Social Survey (ESS) on the Basic
#' Human Values scale (Schwartz, 2007).
#'
#' This dataset includes a circular outcome extracted from the Basic Human
#' Values scale. Note that the extraction of a single circular value as the most
#' salient human value is somewhat of an oversimplification from a theoretical
#' perspective. It is given here because it is nevertheless meaningful as well
#' as useful for illustration purposes.
#'
#' In addition to the circular outcome, some covariates are included. For
#' further details on the variables included, see the ESS documentation.
#'
#' @docType data
#'
#' @usage data(essbhv)
#'
#' @format A data frame with 1690 rows and 10 variables: \describe{
#'   \item{id}{ESS id number}
#'   \item{happy}{Happiness from 0 (extremely unhappy) to 10 (extremely happy).}
#'   \item{rlgdgr}{Self-reported level of how religious one is from 0 (not at all) to 10 (very religious).}
#'   \item{agea}{Age in years.}
#'   \item{edlvenl}{Highest level of education completed on an 18-point scale.}
#'   \item{theta}{Angle on the Basic Human Values scale in radians from -pi to pi.}
#'   \item{thetapos}{theta plus pi. }
#'   \item{thetagrades}{thetapos converted to degrees.}
#'   }
#'
#' @keywords datasets
#'
#' @source \href{https://www.europeansocialsurvey.org/data/}{ESS Data Portal}
#'
#' @references ESS Round 7: European Social Survey Round 7 Data (2014). Data
#'   file edition 2.1. NSD - Norwegian Centre for Research Data, Norway – Data
#'   Archive and distributor of ESS data for ESS ERIC.
#'
#'   Schwartz (2007) Basic human values: theory, methods, and application
"essbhv"

