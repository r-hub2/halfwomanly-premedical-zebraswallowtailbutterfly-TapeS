// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// petterson
double petterson(int sp, double d13);
RcppExport SEXP _TapeS_petterson(SEXP spSEXP, SEXP d13SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type sp(spSEXP);
    Rcpp::traits::input_parameter< double >::type d13(d13SEXP);
    rcpp_result_gen = Rcpp::wrap(petterson(sp, d13));
    return rcpp_result_gen;
END_RCPP
}
// biomass
NumericVector biomass(IntegerVector spp, NumericVector d13, NumericVector d03, NumericVector h);
RcppExport SEXP _TapeS_biomass(SEXP sppSEXP, SEXP d13SEXP, SEXP d03SEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type spp(sppSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type d13(d13SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type d03(d03SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(biomass(spp, d13, d03, h));
    return rcpp_result_gen;
END_RCPP
}
// nsur
NumericMatrix nsur(IntegerVector spp, NumericVector dbh, NumericVector ht, NumericVector sth, NumericVector d03, NumericVector kl);
RcppExport SEXP _TapeS_nsur(SEXP sppSEXP, SEXP dbhSEXP, SEXP htSEXP, SEXP sthSEXP, SEXP d03SEXP, SEXP klSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type spp(sppSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type dbh(dbhSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ht(htSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sth(sthSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type d03(d03SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type kl(klSEXP);
    rcpp_result_gen = Rcpp::wrap(nsur(spp, dbh, ht, sth, d03, kl));
    return rcpp_result_gen;
END_RCPP
}
// nsur2
NumericMatrix nsur2(IntegerVector spp, NumericVector dbh, NumericVector ht);
RcppExport SEXP _TapeS_nsur2(SEXP sppSEXP, SEXP dbhSEXP, SEXP htSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type spp(sppSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type dbh(dbhSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ht(htSEXP);
    rcpp_result_gen = Rcpp::wrap(nsur2(spp, dbh, ht));
    return rcpp_result_gen;
END_RCPP
}
// spline_basis
List spline_basis(NumericVector knots, IntegerVector order, NumericVector xvals, IntegerVector derivs);
RcppExport SEXP _TapeS_spline_basis(SEXP knotsSEXP, SEXP orderSEXP, SEXP xvalsSEXP, SEXP derivsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type knots(knotsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type order(orderSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xvals(xvalsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type derivs(derivsSEXP);
    rcpp_result_gen = Rcpp::wrap(spline_basis(knots, order, xvals, derivs));
    return rcpp_result_gen;
END_RCPP
}
// lmeSKEBLUP
List lmeSKEBLUP(NumericVector xm, NumericVector ym, NumericVector xp, List par, NumericVector RV);
RcppExport SEXP _TapeS_lmeSKEBLUP(SEXP xmSEXP, SEXP ymSEXP, SEXP xpSEXP, SEXP parSEXP, SEXP RVSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type xm(xmSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ym(ymSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xp(xpSEXP);
    Rcpp::traits::input_parameter< List >::type par(parSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type RV(RVSEXP);
    rcpp_result_gen = Rcpp::wrap(lmeSKEBLUP(xm, ym, xp, par, RV));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_TapeS_petterson", (DL_FUNC) &_TapeS_petterson, 2},
    {"_TapeS_biomass", (DL_FUNC) &_TapeS_biomass, 4},
    {"_TapeS_nsur", (DL_FUNC) &_TapeS_nsur, 6},
    {"_TapeS_nsur2", (DL_FUNC) &_TapeS_nsur2, 3},
    {"_TapeS_spline_basis", (DL_FUNC) &_TapeS_spline_basis, 4},
    {"_TapeS_lmeSKEBLUP", (DL_FUNC) &_TapeS_lmeSKEBLUP, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_TapeS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
