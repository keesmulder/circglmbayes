// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/circglmbayes.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// atanLF
vec atanLF(vec x, double r);
RcppExport SEXP _circglmbayes_atanLF(SEXP xSEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(atanLF(x, r));
    return rcpp_result_gen;
END_RCPP
}
// invAtanLF
vec invAtanLF(vec x, double r);
RcppExport SEXP _circglmbayes_invAtanLF(SEXP xSEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(invAtanLF(x, r));
    return rcpp_result_gen;
END_RCPP
}
// rvmc
NumericVector rvmc(int n, double mu, double kp);
RcppExport SEXP _circglmbayes_rvmc(SEXP nSEXP, SEXP muSEXP, SEXP kpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type kp(kpSEXP);
    rcpp_result_gen = Rcpp::wrap(rvmc(n, mu, kp));
    return rcpp_result_gen;
END_RCPP
}
// sampleKappa
vec sampleKappa(double etag, int eta);
RcppExport SEXP _circglmbayes_sampleKappa(SEXP etagSEXP, SEXP etaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type etag(etagSEXP);
    Rcpp::traits::input_parameter< int >::type eta(etaSEXP);
    rcpp_result_gen = Rcpp::wrap(sampleKappa(etag, eta));
    return rcpp_result_gen;
END_RCPP
}
// computeMeanDirection
double computeMeanDirection(vec th);
RcppExport SEXP _circglmbayes_computeMeanDirection(SEXP thSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type th(thSEXP);
    rcpp_result_gen = Rcpp::wrap(computeMeanDirection(th));
    return rcpp_result_gen;
END_RCPP
}
// computeResultantLength
double computeResultantLength(vec th);
RcppExport SEXP _circglmbayes_computeResultantLength(SEXP thSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type th(thSEXP);
    rcpp_result_gen = Rcpp::wrap(computeResultantLength(th));
    return rcpp_result_gen;
END_RCPP
}
// circQuantile
vec circQuantile(arma::vec th, arma::vec q);
RcppExport SEXP _circglmbayes_circQuantile(SEXP thSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type th(thSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(circQuantile(th, q));
    return rcpp_result_gen;
END_RCPP
}
// estimateModeCirc
double estimateModeCirc(NumericVector x, double cip);
RcppExport SEXP _circglmbayes_estimateModeCirc(SEXP xSEXP, SEXP cipSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type cip(cipSEXP);
    rcpp_result_gen = Rcpp::wrap(estimateModeCirc(x, cip));
    return rcpp_result_gen;
END_RCPP
}
// computeHDICirc
NumericVector computeHDICirc(NumericVector x, double cip);
RcppExport SEXP _circglmbayes_computeHDICirc(SEXP xSEXP, SEXP cipSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type cip(cipSEXP);
    rcpp_result_gen = Rcpp::wrap(computeHDICirc(x, cip));
    return rcpp_result_gen;
END_RCPP
}
// estimateMode
double estimateMode(vec x, double cip);
RcppExport SEXP _circglmbayes_estimateMode(SEXP xSEXP, SEXP cipSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type cip(cipSEXP);
    rcpp_result_gen = Rcpp::wrap(estimateMode(x, cip));
    return rcpp_result_gen;
END_RCPP
}
// computeHDI
vec computeHDI(vec x, double cip);
RcppExport SEXP _circglmbayes_computeHDI(SEXP xSEXP, SEXP cipSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type cip(cipSEXP);
    rcpp_result_gen = Rcpp::wrap(computeHDI(x, cip));
    return rcpp_result_gen;
END_RCPP
}
// logProbNormal
vec logProbNormal(vec x, vec mu, vec sd);
RcppExport SEXP _circglmbayes_logProbNormal(SEXP xSEXP, SEXP muSEXP, SEXP sdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< vec >::type sd(sdSEXP);
    rcpp_result_gen = Rcpp::wrap(logProbNormal(x, mu, sd));
    return rcpp_result_gen;
END_RCPP
}
// truncCauchyPdf
double truncCauchyPdf(double x, double m, double w);
RcppExport SEXP _circglmbayes_truncCauchyPdf(SEXP xSEXP, SEXP mSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(truncCauchyPdf(x, m, w));
    return rcpp_result_gen;
END_RCPP
}
// circGLMC
Rcpp::List circGLMC(vec th, mat X, mat D, vec conj_prior, mat bt_prior, vec starting_values, int burnin, int thin, vec bwb, double kappaModeEstBandwith, double CIsize, int Q, double r, bool returnPostSample, int bt_prior_type, bool reparametrize, bool groupMeanComparisons);
RcppExport SEXP _circglmbayes_circGLMC(SEXP thSEXP, SEXP XSEXP, SEXP DSEXP, SEXP conj_priorSEXP, SEXP bt_priorSEXP, SEXP starting_valuesSEXP, SEXP burninSEXP, SEXP thinSEXP, SEXP bwbSEXP, SEXP kappaModeEstBandwithSEXP, SEXP CIsizeSEXP, SEXP QSEXP, SEXP rSEXP, SEXP returnPostSampleSEXP, SEXP bt_prior_typeSEXP, SEXP reparametrizeSEXP, SEXP groupMeanComparisonsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vec >::type th(thSEXP);
    Rcpp::traits::input_parameter< mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< mat >::type D(DSEXP);
    Rcpp::traits::input_parameter< vec >::type conj_prior(conj_priorSEXP);
    Rcpp::traits::input_parameter< mat >::type bt_prior(bt_priorSEXP);
    Rcpp::traits::input_parameter< vec >::type starting_values(starting_valuesSEXP);
    Rcpp::traits::input_parameter< int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< int >::type thin(thinSEXP);
    Rcpp::traits::input_parameter< vec >::type bwb(bwbSEXP);
    Rcpp::traits::input_parameter< double >::type kappaModeEstBandwith(kappaModeEstBandwithSEXP);
    Rcpp::traits::input_parameter< double >::type CIsize(CIsizeSEXP);
    Rcpp::traits::input_parameter< int >::type Q(QSEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    Rcpp::traits::input_parameter< bool >::type returnPostSample(returnPostSampleSEXP);
    Rcpp::traits::input_parameter< int >::type bt_prior_type(bt_prior_typeSEXP);
    Rcpp::traits::input_parameter< bool >::type reparametrize(reparametrizeSEXP);
    Rcpp::traits::input_parameter< bool >::type groupMeanComparisons(groupMeanComparisonsSEXP);
    rcpp_result_gen = Rcpp::wrap(circGLMC(th, X, D, conj_prior, bt_prior, starting_values, burnin, thin, bwb, kappaModeEstBandwith, CIsize, Q, r, returnPostSample, bt_prior_type, reparametrize, groupMeanComparisons));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_circglmbayes_atanLF", (DL_FUNC) &_circglmbayes_atanLF, 2},
    {"_circglmbayes_invAtanLF", (DL_FUNC) &_circglmbayes_invAtanLF, 2},
    {"_circglmbayes_rvmc", (DL_FUNC) &_circglmbayes_rvmc, 3},
    {"_circglmbayes_sampleKappa", (DL_FUNC) &_circglmbayes_sampleKappa, 2},
    {"_circglmbayes_computeMeanDirection", (DL_FUNC) &_circglmbayes_computeMeanDirection, 1},
    {"_circglmbayes_computeResultantLength", (DL_FUNC) &_circglmbayes_computeResultantLength, 1},
    {"_circglmbayes_circQuantile", (DL_FUNC) &_circglmbayes_circQuantile, 2},
    {"_circglmbayes_estimateModeCirc", (DL_FUNC) &_circglmbayes_estimateModeCirc, 2},
    {"_circglmbayes_computeHDICirc", (DL_FUNC) &_circglmbayes_computeHDICirc, 2},
    {"_circglmbayes_estimateMode", (DL_FUNC) &_circglmbayes_estimateMode, 2},
    {"_circglmbayes_computeHDI", (DL_FUNC) &_circglmbayes_computeHDI, 2},
    {"_circglmbayes_logProbNormal", (DL_FUNC) &_circglmbayes_logProbNormal, 3},
    {"_circglmbayes_truncCauchyPdf", (DL_FUNC) &_circglmbayes_truncCauchyPdf, 3},
    {"_circglmbayes_circGLMC", (DL_FUNC) &_circglmbayes_circGLMC, 17},
    {NULL, NULL, 0}
};

RcppExport void R_init_circglmbayes(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
