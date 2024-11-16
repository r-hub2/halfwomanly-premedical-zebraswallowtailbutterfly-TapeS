#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

//' height estimation
//'
//' height estimation based on diameter in breast height and species using a
//' Petterson-function
//'
//' @param sp vector of species code for biomass function from interval [1;18];
//' see \code{\link{BaMap}} for mapping of species model codes
//' @param d13 vector of diameter in breast height; in centimeter
//' @return a scalar: tree height
//' @export
// [[Rcpp::export]]
double petterson(int sp, double d13){
  double a[18] = {0.274072, 0.277236, 0.259403, 0.297216, 0.276316,
                  0.293966, 0.315672, 0.300388, 0.332180, 0.314044,
                  0.317395, 0.311777, 0.322445, 0.326550, 0.321719,
                  0.322696, 0.280639, 0.345238};
  double b[18] = {2.220311, 2.397463, 2.928833, 1.986882, 2.459936,
                  1.768945, 1.633350, 1.535308, 1.281238, 1.532003,
                  1.344200, 1.729393, 1.490083, 1.292453, 1.575659,
                  1.496584, 2.402878, 1.695635};
  return (1.3 + 1 / pow(a[sp-1] + b[sp-1] / d13, 3.0));
}

//' Prediction of above-ground biomass according to NFI-functions
//'
//' Prediction of total above-ground biomass for trees defined via species, dbh,
//' d03 and height
//' @param spp vector of species code for biomass function [1;18]
//' @param d13 vector of diameter in breast height in centimeter
//' @param d03 vector of diameter in 30\% of tree height in centimeter
//' @param h vector of height of trees
//' @return a vector of total above-ground biomass
//' @details code taken from BDAT (Koeff.f).
//' @references
//' Riedel, T. and G. Kaendler (2017). "Nationale Treibhausgasberichterstattung:
//' Neue Funktionen zur Schätzung der oberirdischen Biomasse am Einzelbaum."
//' Forstarchiv 88(2): 31-38.
//' @export
// [[Rcpp::export]]
NumericVector biomass(IntegerVector spp, NumericVector d13, NumericVector d03,
                      NumericVector h) {


  int nx = spp.length();  // number of elements to process
  int bax;                // index of species code to be used to access coefs
                          // in fortran code this is 'ba_nr'
  double d03_os;          // obere schranke für d03
  double h_os;            // obere schranke für h
  NumericVector bm(nx);   // resulting vector bm=biomass

  // coefficients
  double b0_h[18] = {0.230589, 0.230589, 0.230589, 0.230589, 0.230589,
                        0.049402, 0.049402, 0.049402, 0.049402, 0.049402,
                        0.049402, 0.049402, 0.049402, 0.049402, 0.049402,
                        0.049402, 0.049402, 0.049402};
  double b1_h[18] = {2.201010, 2.201010, 2.201010, 2.201010, 2.201010,
                        2.549462, 2.549462, 2.549462, 2.549462, 2.549462,
                        2.549462, 2.549462, 2.549462, 2.549462, 2.549462,
                        2.549462, 2.549462, 2.549462};
  double Bo[18] = {0.410799, 0.410799, 0.410799, 0.410799, 0.410799,
                      0.096436, 0.096436, 0.096436, 0.096436, 0.096436,
                      0.096436, 0.096436, 0.096436, 0.096436, 0.096436,
                      0.096436, 0.096436, 0.096436};
  double B_us[18] = {26.631220, 19.122312, 19.117109, 19.999435, 28.339055,
                        33.223278, 28.947819, 33.748928, 38.949475, 30.795513,
                        33.090857, 22.414603, 41.656504, 33.975677, 28.289641,
                        19.888007, 16.861011, 21.785135};
  double b3_poly[18] = {0.013696, -0.000710, -0.002753, 0.009158, 0.027686,
                           0.011621, 0.015009, 0.026795, 0.023861, 0.024696,
                           0.029672, 0.006530, 0.041933, 0.025661, 0.019265,
                           0.001053, -0.005509, 0.005165};
  double b0[18] = {0.752848, 0.063943, 0.054267, 0.337783, 0.141263,
                      0.167871, 0.094279, 0.140114, 0.348054, 0.207830,
                      0.206695, 0.110756, 0.627748, 0.191710, 0.210419,
                      0.138305, 0.272778, 0.299783};
  double b1[18] = {2.849849, 4.112771, 5.531211, 2.840552, 4.419252,
                      6.254522, 10.269984, 4.917100, 5.512190, 4.379202,
                      3.944298, 5.176043,	3.971093, 4.202065, 4.387211,
                      4.463595, 4.192402, 4.755305};
  double b2[18] = {6.030355, 6.692779, 9.259898, 6.349639, 6.294424,
                      6.647523, 8.138936, 6.285623, 5.829990, 6.294122,
                      6.113504, 6.965050, 5.135899, 6.656169, 5.439829,
                      5.947719, 5.962976, 6.365869};
  double b3[18] = {0.621878, 1.058685, 0.890083, 0.627555, 0.883409,
                      0.807451, 0.558449, 0.905214, 0.884464, 0.859312,
                      0.862962, 0.852889, 0.870941, 0.779942, 0.806097,
                      0.850972, 0.810312, 0.786933};
  double k1[18] = {42.00, 6.00,   5.50, 18.00, 6.40, 11.00,
                      400.00, 9.00, 157.60, 10.00, 9.40, 8.70,
                      67.60, 8.20, 68.40, 11.10, 13.70, 13.20};
  double k2[18] = {24.0, 62.70, 139.20, 23.00, 68.60, 135.0,
                      8.0, 70.60, 13.00, 61.50, 46.30, 90.90,
                      18.60, 53.10, 10.10, 50.70, 66.80, 85.80};
  double d13_os[18] = {69.00, 90.00, 82.00, 59.00, 70.00, 86.0,
                          94.00, 77.00, 55.00, 68.00, 69.00, 84.0,
                          75.00, 85.00, 53.00, 58.00, 113.00, 102.0};
  double c0[18] = {1.078434, 1.173719, 1.068712, 0.890090, 1.200444,
                      0.840135, 0.876335, 0.828197, 0.801798, 0.864267,
                      0.748943, 0.830802, 0.980605, 1.035449, 1.030980,
                      0.936909, 0.867197, 0.828710};
  double c1[18] = {0.912040, 0.898113, 0.906065, 0.957474, 0.880278,
                      0.989704, 0.982788, 0.991109, 1.000543, 0.975139,
                      1.021867, 0.990512, 0.935437, 0.926549, 0.904414,
                      0.944473, 0.961536, 0.986861};

  // processing each observation
  for(int i = 0; i < nx; i++) {

    bax = spp[i] - 1;      // vector index from 0 to 17 = 18 elements
    if(bax > 17) bax = 1;   // assure the range of validity

    d03_os = d03[i] + c0[bax] * pow(d13_os[bax], c1[bax]) - c0[bax] *
      pow(d13[i], c1[bax]);

    h_os = h[i] + petterson(spp[i], d13_os[bax]) - petterson(spp[i], d13[i]);

    if(h[i] < 1.3){

      bm[i] = b0_h[bax] * pow(h[i], b1_h[bax]);

    } else if(d13[i] < 10) {

      bm[i] = Bo[bax] + ((B_us[bax] - Bo[bax]) / (10*10) +
        b3_poly[bax] * (d13[i] - 10)) * (d13[i] * d13[i]);

    } else if(d13[i] < d13_os[bax]) {

      bm[i] = b0[bax] * exp( b1[bax]  * d13[i] / (d13[i] + k1[bax]) ) *
        exp( b2[bax] * d03[i] / (d03[i] + k2[bax]) ) * pow(h[i], b3[bax]);

    } else {

      bm[i] = b0[bax] * exp( b1[bax] * d13_os[bax] / (d13_os[bax] + k1[bax]) ) *
        exp( b2[bax] * d03_os / (d03_os + k2[bax]) ) * pow(h_os, b3[bax]);

      bm[i] = bm[i] * (1 + b1[bax] * k1[bax] /
        pow(d13_os[bax] + k1[bax], 2) * (d13[i] - d13_os[bax]) +
        b2[bax] * k2[bax] / pow(d03_os + k2[bax], 2) * (d03[i] - d03_os) +
        b3[bax] / h_os * (h[i] - h_os));

    }
  }

  return bm;
}

//' Component biomass functions
//'
//' evaluation of the component biomass functions fit by nonlinear seemingly
//' unrelated regression (NSUR) to estimate absolute or relative component mass
//'
//' @param spp vector of species code for biomass component function of interval
//' [1;8]; see \code{\link{BaMap}} for mapping of species model codes
//' @param dbh vector of diameter in breast height; in centimeter
//' @param ht vector of tree heights, in meter
//' @param sth vector of stump heights, in meter
//' @param d03 vector if diameter in 30\% of tree height, in centimeter
//' @param kl vector of crown length, i.e. tree height minus height of crown base, in meter
//' @return a numeric matrix holding component biomass
//' @details function to calculate component biomass; functions fitted using
//' same methodology as in Vonderach et al. (2018) with slightly updated
//' parameters as in Vonderach and Kändler (2021); species mapping as in
//' \code{TapeS::BaMap(, type=7)};
//' @references Vonderach, C., G. Kändler and C. F. Dormann (2018).
//' "Consistent set of additive biomass functions for eight tree species in
//' Germany fit by nonlinear seemingly unrelated regression."
//' Annals of Forest Science 75(2): 49.
//' \doi{10.1007/s13595-018-0728-4}
//'
//' Vonderach, C. and G. Kändler (2021). Neuentwicklung von Schaftkurven- und
//' Biomassemodellen für die Bundeswaldinventur auf Basis des TapeR-Pakets -
//' Abschlussbericht zum Projekt BWI-TapeR. Freiburg: 150p.
//' @examples
//' nsur(spp = c(1, 6),
//'      dbh = c(30, 30),
//'      ht = c(25, 27),
//'      sth = c(0.25, 0.27),
//'      d03 = c(27, 27),
//'      kl = .7*c(25, 27))
//' @export
// [[Rcpp::export]]
NumericMatrix nsur(IntegerVector spp, NumericVector dbh, NumericVector ht,
                   NumericVector sth, NumericVector d03, NumericVector kl){
  int n = spp.size();
  //check input vector spp
  for(int i = 0; i < n; ++i) {

    if((spp[i] > 8 || spp[i] < 1 || IntegerVector::is_na(spp[i]))){
      spp[i] = 1;
    }

  }

  // defining parameters for biomass equation
  // a + b * BHD^c * Hoehe^d * D03^e * Stockhoehe^f * KL^g * DH^h
  // parameter array each with 8 species (rows) and 6 components (cols)
  // rows=(fi, ta, kie, dgl, bu, ei, bah, es)
  // cols=(stump, stump bark, coarse wood, coarse wood bark, small wood, needles)
  double a[8][6] = {
    {0, 0, 0.6999, 0, 0, -1.9558},
    {0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 5.3977, 0},
    {0, 0, -4.8678, 0, 0, 0},
    {0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 7.2819, 0},
    {0, 0, 0, -6.9386, 0, 0}
  };
  double b[8][6] = {
    {0.0211, 0.0058, 0.0125, 0.0034, 0.0435, 0.2496},
    {0.0151, 0.0039, 0.0062, 0.0028, 0.0408, 0.1621},
    {0.0623, 0.0075, 0.0192, 0.0042, 0.265, 0.2206},
    {0.0174, 0.0015, 0.0102, 0.0012, 0.0081, 0.0427},
    {0.0303, 0.0039, 0.0214, 0.0016, 0.4248, 0},
    {0.0359, 0.0153, 0.0156, 0.0052, 0.3631, 0},
    {0.0487, 0.0095, 0.0178, 0.0045, 1e-04, 0},
    {0.0126, 0.4908, 0.0146, 1, 0.1298, 0}
  };
  double c[8][6] = {
    {2.1363, 1.7085, 0.8808, 0.4272, 1.5836, 1.2606},
    {2.1968, 2.0958, 0.4735, 0.3706, 2.5348, 1.3889},
    {1.9326, 1.8205, 1.3807, 0.7706, 0, 0},
    {2.2085, 1.8333, 0.4985, 0, 0, 0},
    {2.1528, 1.9243, 1.4721, 1.2836, 1.5779, 0},
    {2.0709, 1.7996, 1.3635, 1.2974, 1.9356, 0},
    {2.0455, 1.8586, 1.5682, 2.0598, 0, 0},
    {2.4253, 2.3043, 2.0281, 0, 1.5338, 0}
  };
  double d[8][6] = {
    {0, 0, 1.2265, 1.0682, -0.7975, -0.8523},
    {0, 0, 1.7171, 1.1918, -0.6412, 0.1598},
    {0, 0, 0.9548, 0.7427, -0.9237, -1.4946},
    {0, 0.4247, 1.2845, 1.409, -1.7723, -1.3532},
    {0, 0, 0.9363, 1.0719, -0.9679, 0},
    {0, 0, 0.9419, 0.9111, -0.8488, 0},
    {0, 0, 0.8855, 0.6794, 0, 0},
    {0, -1.7406, 1.0797, 0, 0, 0}
  };
  double e[8][6] = {
    {0, 0.2019, 0.9769, 1.2822, 1.2901, 0.9867},
    {0, 0, 1.1382, 1.4819, 0, 0},
    {0, 0, 0.6078, 1.2371, 2.0524, 2.3028},
    {0, 0, 1.4067, 1.7932, 3.1899, 2.6678},
    {0, 0, 0.6291, 0.689, 0.5343, 0},
    {0, 0, 0.805, 0.6662, 0, 0},
    {0, 0, 0.5722, 0, 0, 0},
    {0, 0, -0.0027, 3.2731, 0, 0}
  };
  double f[8][6] = {
    {0.6076, 0.8383, 0, 0, 0, 0},
    {0.761, 0.7683, 0, 0, 0, 0},
    {1.0437, 0.8745, 0, 0, 0, 0},
    {0.7966, 0.7928, 0, 0, 0, 0},
    {0.7959, 0.8059, 0, 0, 0, 0},
    {0.7803, 0.9041, 0, 0, 0, 0},
    {0.814, 0.8242, 0, 0, 0, 0},
    {0.7427, 0.7169, 0, 0, 0, 0}
  };
  double g[8][6] = {
    {0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0.3444, 0.3374},
    {0, 0, 0, 0, 0.5434, 0.5892},
    {0, 0, 0, 0, 1.4282, 0.7153},
    {0, 0, 0, 0, 0.543, 0},
    {0, 0, 0, 0, 0.5139, 0},
    {0, 0, 0, 0, 2.3559, 0},
    {0, 0, 0, 0, 0.4926, 0}
  };
  double h[8][6] = {
    {0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 1.0556, 0},
    {0, 0, 0, -0.9729, 0, 0}
  };

  int nx = spp.length();  // number of elements to process
  IntegerVector idx = seq_along(spp);
  NumericVector DH = dbh * (ht - sth); //only for species 'es'
  NumericMatrix res(spp.length(), 7);


  res(_, 0) = idx; // referencing given observation for later purpose
  for(int i=0; i<nx; i++){
    for(int j=0; j<6; j++){
      // Rcpp::Rcout << a[ spp[i]-1 ][ j ] << std::endl;
      res(i, j+1) = a[ spp[i]-1 ][ j ] +
        b[ spp[i]-1 ][ j ] *
        pow(dbh[i], c[ spp[i]-1 ][ j ]) *
        pow((ht[i] - sth[i]), d[ spp[i]-1 ][ j ]) *
        pow(d03[i], e[ spp[i]-1 ][ j ]) *
        pow(sth[i], f[ spp[i]-1 ][ j ]) *
        pow(kl[i], g[ spp[i]-1 ][ j ]) *
        pow(DH[i], h[ spp[i]-1 ][ j ]);
    }
  }
  // adding column names
  colnames(res) = CharacterVector::create("id", "stw", "stb", "sw", "sb",
           "fwb", "ndl");

  return res;
}

//' Component biomass functions
//'
//' evaluation of the component biomass functions fit by nonlinear seemingly
//' unrelated regression (NSUR) to estimate absolute or relative component mass
//'
//' @param spp vector of species code for biomass component function of interval
//' [1;8]; see \code{\link{BaMap}} for mapping of species model codes
//' @param dbh vector of diameter in breast height; in centimeter
//' @param ht vector of tree heights, in meter
//' @return a numeric matrix holding component biomass
//' @details simple function from Vonderach et al. (2018) to calculate component
//' biomass; species mapping as in \code{TapeS::BaMap(, type=7)}
//' @references Vonderach, C., G. Kändler and C. F. Dormann (2018).
//' "Consistent set of additive biomass functions for eight tree species in
//' Germany fit by nonlinear seemingly unrelated regression."
//' Annals of Forest Science 75(2): 49.
//' \doi{10.1007/s13595-018-0728-4}
//'
//' @examples
//' nsur2(spp = c(1, 6),
//'       dbh = c(30, 30),
//'       ht = c(25, 27))
//' @export
// [[Rcpp::export]]
NumericMatrix nsur2(IntegerVector spp, NumericVector dbh, NumericVector ht){
  int n = spp.size();
  //check input vector spp
  for(int i = 0; i < n; ++i) {

    if((spp[i] > 8 || spp[i] < 1 || IntegerVector::is_na(spp[i]))){
      spp[i] = 1;
    }

  }

  // defining parameters for biomass equation
  // a + b * BHD^c * Hoehe^d * Stockhoehe^e
  // parameter array each with 8 species (rows) and 6 components (cols)
  // rows=(fi, ta, kie, dgl, bu, ei, bah, es)
  // cols=(stump, stump bark, coarse wood, coarse wood bark, small wood, needles)
  double a[8][6] = { // intcpt
    {0, -0.0128, 0, 0, 1.8472, -1.6847},  //fi
    {0, 0, 0, 0, 0, 0},                   //ta
    {0, 0, 0, 0, 0, 0},                   //kie
    {0, 0, 0, 0, 0, -1.8821},             //dgl
    {0, 0, -4.6332, 0, 0, 0},             //bu
    {0, 0, -3.9731, 0, 0, 0},             //ei
    {0, 0, 0, 0, 0, 0},                   //bah
    {0, 0, 0, 0, 0, 0}                    //es
  };
  double b[8][6] = { //scale
    {0.0220, 0.0067, 0.0142, 0.0038, 0.0243, 0.2850}, //fi
    {0.0121, 0.0036, 0.0046, 0.0019, 0.0273, 0.1071}, //ta
    {0.0624, 0.0077, 0.0173, 0.0055, 0.1316, 0.1484}, //kie
    {0.0186, 0.0032, 0.0131, 0.0018, 0.2784, 0.2749}, //dgl
    {0.0315, 0.0040, 0.0190, 0.0013, 0.4006, 0},      //bu
    {0.0363, 0.0158, 0.0223, 0.0043, 0.3028, 0},      //ei
    {0.0596, 0.0121, 0.0170, 0.0040, 0.0426, 0},      //bah
    {0.0111, 0.0367, 0.0220, 0.0006, 0.0966, 0}       //es
  };
  double c[8][6] = { // dbh
    {2.1212, 1.7268, 1.7414, 1.6076, 2.9671, 2.1173}, //fi
    {2.2645, 2.1225, 1.1917, 1.4458, 2.2573, 1.6952}, //ta
    {1.9322, 1.8127, 2.0072, 2.0108, 2.4440, 2.3320}, //kie
    {2.1850, 2.0357, 1.9299, 1.9099, 3.1276, 2.4833}, //dgl
    {2.1447, 1.9184, 2.0861, 1.9596, 2.3211, 0},      //bu
    {2.0657, 1.7910, 2.1012, 1.9437, 2.1945, 0},      //ei
    {1.9934, 1.8062, 2.0710, 2.0287, 2.1166, 0},      //bah
    {2.4593, 1.4515, 2.0683, 1.7301, 1.9989, 0}       //es
  };
  double d[8][6] = { // ht
    {0, 0, 1.2401, 1.0528, -0.8183, -0.8334}, //fi
    {0, 0, 2.2103, 1.6780, 0, 0},             //ta
    {0, 0, 0.9140, 0.5374, -0.9490, -1.2026}, //kie
    {0, 0, 1.0715, 1.0306, -1.7984, -1.3051}, //dgl
    {0, 0, 0.9502, 1.0990, -0.7636, 0},       //bu
    {0, 0, 0.8647, 0.9576, -0.6882, 0},       //ei
    {0, 0, 0.9317, 0.7446, 0, 0},             //bah
    {0, 0, 0.9050, 1.7301, 0, 0}              //es
  };
  double e[8][6] = { //stump-ht
    {0.6056, 0.5947, 0, 0, 0, 0}, //fi
    {0.7596, 0.7856, 0, 0, 0, 0}, //ta
    {1.0414, 0.8732, 0, 0, 0, 0}, //kie
    {0.7723, 0.7621, 0, 0, 0, 0}, //dgl
    {0.7980, 0.8076, 0, 0, 0, 0}, //bu
    {0.7721, 0.9032, 0, 0, 0, 0}, //ei
    {0.8314, 0.8178, 0, 0, 0, 0}, //bah
    {0.7579, 0.7298, 0, 0, 0, 0}  //es
  };


  int nx = spp.length();  // number of elements to process
  IntegerVector idx = seq_along(spp);
  NumericVector sth = 0.007 * ht; // according to the fitting data see devel/check-sth-kro.r
  NumericMatrix res(spp.length(), 7);


  res(_, 0) = idx; // referencing given observation for later purpose
  for(int i=0; i<nx; i++){
    for(int j=0; j<6; j++){
      // Rcpp::Rcout << a[ spp[i]-1 ][ j ] << std::endl;
      res(i, j+1) = a[ spp[i]-1 ][ j ] +
        b[ spp[i]-1 ][ j ] *
        pow(dbh[i], c[ spp[i]-1 ][ j ]) *
        pow(ht[i], d[ spp[i]-1 ][ j ]) *
        pow(sth[i], e[ spp[i]-1 ][ j ]);
    }
  }
  // adding column names
  colnames(res) = CharacterVector::create("id", "stw", "stb", "sw", "sb",
           "fwb", "ndl");

  return res;
}
