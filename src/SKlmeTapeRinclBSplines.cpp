/* This code originates from R-package splines, i.e. splines.c
 * and was slightly modified to be used in c-context using Rcpp
 * Modifier: christian.vonderach@forst.bwl.de
 * Date: 23.06.2020
 * All license stuff still applies.
 */

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

typedef struct spl_struct {
  int order,			/* order of the spline */
    ordm1,			/* order - 1 (3 for cubic splines) */
    nknots,			/* number of knots */
    curs,			/* current position in knots vector */
    boundary;		/* must have knots[curs] <= x < knots[curs+1] */
    /* except for the boundary case */

    double *ldel,		/* differences from knots on the left */
    *rdel,			/* differences from knots on the right */
    *knots,			/* knot vector */
    *coeff,			/* coefficients */
    *a;			/* scratch array */
} *splPTR;

/* set sp->curs to the index of the first knot position > x.
Special handling for x == sp->knots[sp->nknots - sp-order + 1] */
  static int
set_cursor(splPTR sp, double x)
{
  int i;
  /* don't assume x's are sorted */

    sp->curs = -1; /* Wall */
    sp->boundary = 0;
  for (i = 0; i < sp->nknots; i++) {
    if (sp->knots[i] >= x) sp->curs = i;
    if (sp->knots[i] > x) break;
  }
  if (sp->curs > sp->nknots - sp->order) {
    int lastLegit = sp->nknots - sp->order;
    if (x == sp->knots[lastLegit]) {
      sp->boundary = 1; sp->curs = lastLegit;
    }
  }
  return sp->curs;
}

  static void
diff_table(splPTR sp, double x, int ndiff)
{
  int i;
  for (i = 0; i < ndiff; i++) {
    sp->rdel[i] = sp->knots[sp->curs + i] - x;
    sp->ldel[i] = x - sp->knots[sp->curs - (i + 1)];
  }
}


/* fast evaluation of basis functions */
  static void
basis_funcs(splPTR sp, double x, double *b)
{
  diff_table(sp, x, sp->ordm1);
  b[0] = 1.;
  for (int j = 1; j <= sp->ordm1; j++) {
    double saved = 0.;
    for (int r = 0; r < j; r++) { // do not divide by zero
      double den = sp->rdel[r] + sp->ldel[j - 1 - r];
      if(den != 0) {
        double term = b[r]/den;
        b[r] = saved + sp->rdel[r] * term;
        saved = sp->ldel[j - 1 - r] * term;
      } else {
        if(r != 0 || sp->rdel[r] != 0.)
          b[r] = saved;
        saved = 0.;
      }
    }
    b[j] = saved;
  }
}

/* "slow" evaluation of (derivative of) basis functions */
static double
  evaluate(splPTR sp, double x, int nder)
  {
    double *lpt, *rpt, *apt, *ti = sp->knots + sp->curs;
    int inner, outer = sp->ordm1;

    if (sp->boundary && nder == sp->ordm1) { /* value is arbitrary */
return 0.0;
    }
    while(nder--) {  // FIXME: divides by zero
      for(inner = outer, apt = sp->a, lpt = ti - outer; inner--; apt++, lpt++)
        *apt = outer * (*(apt + 1) - *apt)/(*(lpt + outer) - *lpt);
      outer--;
    }
    diff_table(sp, x, outer);
    while(outer--)
      for(apt = sp->a, lpt = sp->ldel + outer, rpt = sp->rdel, inner = outer + 1;
          inner--; lpt--, rpt++, apt++)
        // FIXME: divides by zero
        *apt = (*(apt + 1) * *lpt + *apt * *rpt)/(*rpt + *lpt);
    return sp->a[0];
  }

// [[Rcpp::export]]
  List
spline_basis(NumericVector knots, IntegerVector order, NumericVector xvals, IntegerVector derivs)
{
  /* evaluate the non-zero B-spline basis functions (or their derivatives)
  * at xvals.  */

  double *kk = REAL(knots);
  int nk = knots.length();
  int ord = as<int>(order);
  double *xx = REAL(xvals);
  int nx = xvals.length();
  int *ders = INTEGER(derivs);
  int nd = derivs.length();

  splPTR sp = (struct spl_struct *) R_alloc(1, sizeof(struct spl_struct));
  /* fill sp : */
  sp->order = ord;
  sp->ordm1 = ord - 1;
  sp->rdel = (double *) R_alloc(sp->ordm1, sizeof(double));
  sp->ldel = (double *) R_alloc(sp->ordm1, sizeof(double));
  sp->knots = kk;
  sp->nknots = nk;
  sp->a = (double *) R_alloc(ord, sizeof(double));
  NumericMatrix val (ord, nx);
  IntegerVector offsets(nx);
  double *valM = REAL(val);
  int *ioff = INTEGER(offsets);

  for(int i = 0; i < nx; i++) {
    set_cursor(sp, xx[i]);
    // ==> io  \in {0,..,nk} is the knot-interval "number"
    int io = ioff[i] = sp->curs - ord;
    int der_i = ders[i % nd];
    if (io < 0 || io > nk) {
      for (int j = 0; j < ord; j++) {
        valM[i * ord + j] = R_NaN;
      }
    } else if (der_i > 0) { /* slow method for derivatives */
      if (der_i >= ord) {
        if(nd == 1) {
          stop("derivs should be in {0,..,ord-1}");
          // error(_("derivs = %d >= ord = %d, but should be in {0,..,ord-1}"),
          //       der_i, ord);
        } else {
          stop("derivs should be in {0,..,ord-1}");
          // error(_("derivs[%d] = %d >= ord = %d, but should be in {0,..,ord-1}"),
          //       i+1, der_i, ord);
        }
      }
      for(int ii = 0; ii < ord; ii++) {
        for(int j = 0; j < ord; j++) sp->a[j] = 0;
        sp->a[ii] = 1;
        valM[i * ord + ii] =
          evaluate(sp, xx[i], der_i);
      }
    } else {		/* fast method for value */
        basis_funcs(sp, xx[i], valM + i * ord);
    }
  }
  /* transform result into DesignMatrix for further processing
   * TODO: case of outer.ok = TRUE not yet implemented, but should not occur!
   * check that knots zero and one are included in knot vector
   */
  int ncoef = nk - ord;
  NumericMatrix design (nx, ncoef); /* initialised with 0 by default as desired */
  for(int i=0; i < nx; i++){
    for(int j=0; j < ord; j++){
      design(i, j + offsets(i) ) = val(j, i);
    }
  }
  NumericMatrix::Sub sub = design( _ , Range(0, ncoef - 2) );

  return List::create(Named("val") = val
                      , Named("Offset") = offsets
                      , Named("design") = design
                      , Named("ncoef") = ncoef
                      , Named("subDesign") = sub
                      );
}


/* diameter prediction E[d] for TapeR-object
 * using BSpline Matrix from above
 */
//' diameter prediction E[d] for TapeR-object
//'
//' Prediction diameter (no variances) for given tree and TapeR-object using
//' BSpline Matrix all in C++
//' @param xm relative height of measured diameter
//' @param ym measured diameter for calibration
//' @param xp relative height for which diameter prediction is required
//' @param par a TapeR-object (including padded knots vector)
//' @param RV numeric vector holding assumed residual variance for each observation
//' @return a list holding several elements, perspectively only the estimated diameter
//' @details code implementation in C++ following the code base of TapeR. Bspline
//'   matrix code taken from R-package splines to avoid the need of calling R from C.
//' @export
// [[Rcpp::export]]
  List
lmeSKEBLUP(NumericVector xm, NumericVector ym, NumericVector xp, List par, NumericVector RV){

  NumericVector x_k = xm;
  arma::vec y_k = Rcpp::as<vec>(ym);
  arma::vec resvar = Rcpp::as<vec>(RV);
  arma::vec b_fix = Rcpp::as<vec>(par["b_fix"]);
  // double sig2_eps = Rcpp::as<double>(par["sig2_eps"]);
  NumericVector dfRes = par["dfRes"];
  NumericMatrix KOVb_fix = par["KOVb_fix"];
  List BSm;
  List BSp;

  NumericVector a = as<NumericVector>(par["pad_knt_x"]);
  IntegerVector b = as<IntegerVector>(par["ord_x"]);
  NumericVector c = xm; // now take data to predict random effect
  IntegerVector d(1); // derivative taken to be zero
  BSm = spline_basis(a, b, c, d);
  mat BSX = as<mat>(BSm["subDesign"]);

  NumericVector az = as<NumericVector>(par["pad_knt_z"]);
  IntegerVector bz = as<IntegerVector>(par["ord_z"]);
  BSm = spline_basis(az, bz, c, d);
  mat BSZ = as<mat>(BSm["subDesign"]);

  // first arma conversion
  mat KOVb = as<mat>(par["KOVb_rnd"]);
  arma::mat Z_KOVb_Zt_k = BSZ * KOVb * trans(BSZ);

  // build R_k
  arma::mat R_k = diagmat(resvar);

  // extend Z_KOVb_Zt_k
  arma::mat KOVinv_y_k = arma::inv(Z_KOVb_Zt_k + R_k);

  // estimate random effects
  arma::vec EBLUP_b_k = KOVb * arma::trans(BSZ) * KOVinv_y_k * (y_k - BSX * b_fix);

  // estimate diameter/s
  // BSp = f(Named("x")=xp, Named("par.lme")=par);
  c = xp; // now take data to predict diameter
  BSp = spline_basis(a, b, c, d);
  mat BSpX = as<mat>(BSp["subDesign"]);
  BSp = spline_basis(az, bz, c, d);
  mat BSpZ = as<mat>(BSp["subDesign"]);
  arma::mat EBLUP_y_kh = BSpX * b_fix + BSpZ * EBLUP_b_k;

  return List::create(Named("yp") = EBLUP_y_kh
  //                       , Named("b_fix") = b_fix
  //                       , Named("b_rnd") = EBLUP_b_k
  //                       , Named("a") = a
  //                       , Named("az") = az
  //                       , Named("b") = b
  //                       , Named("bz") = bz
  //                       , Named("c") = c
  //                       , Named("d") = d
  //                       , Named("BSX") = BSX
  //                       , Named("BSZ") = BSZ
  //                       , Named("BSpX") = BSpX
  //                       , Named("BSpZ") = BSpZ
  //                       , Named("sig2") = sig2_eps
  //                       , Named("R_k") = R_k
  //                       , Named("KOVb") = KOVb
  //                       , Named("Z_KOVb_Zt_k") = Z_KOVb_Zt_k
  //                       , Named("KOVinv_y_k") = KOVinv_y_k
  );
}
