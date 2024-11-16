#' @title Tariff for taper form
#' @description evaluates tariff functions to estimate taper form, i.e. quotient
#' of d03 by d005
#' @param spp species code of \code{tprSpeciesCode}
#' @param Dbh diameter of considered tree at 1.3m above ground [cm]
#' @param Ht tree height of considered tree [m]
#' @param inv indicator for inventory (0=TapeS taper curve models,
#' 1=NFI1, 2=NsoG, 3=IS08, 4=NFI3, 5=BDAT)
#' @return quotient of d03 / d005 [unitless]
#' @references c.f. rBDAT::getForm respectively BDAT source code FormTarif.f
#' @import TapeR
## #' @importFrom rBDAT getForm
#' @export
#' @examples
#' ## dont't run
#' spp <- 15
#' Dbh <- 30
#' Ht <- 27
#' FormTariff(spp, Dbh, Ht, 0)
#' FormTariff(spp=c(1:2), Dbh=c(30, 30), Ht=c(27, 24), inv=0)
#' if(require("rBDAT")){
#'   FormTariff(spp, Dbh, Ht, 0)
#'   rBDAT::getForm(list(spp=spp, D1=Dbh, H1=1.3, H=Ht), inv=0) # different taper curves!
#'   FormTariff(spp, Dbh, Ht, 1)
#'   rBDAT::getForm(list(spp=spp, D1=Dbh, H1=1.3, H=Ht), inv=1) # identical
#'   FormTariff(spp, Dbh, Ht, 2)
#'   rBDAT::getForm(list(spp=spp, D1=Dbh, H1=1.3, H=Ht), inv=2) # identical
#'   FormTariff(spp, Dbh, Ht, 3)
#'   rBDAT::getForm(list(spp=spp, D1=Dbh, H1=1.3, H=Ht), inv=3) # identical
#'   FormTariff(spp, Dbh, Ht, 4)
#'   rBDAT::getForm(list(spp=spp, D1=Dbh, H1=1.3, H=Ht), inv=4) # identical
#' }
#'


FormTariff <- function(spp, Dbh, Ht, inv){
  EQP <- getEQP()
  EQ03uG <- 0.4
  EQ03oG <- 0.98
  phi <- 0.001
  Tarif <- as.integer(inv)

  if(Tarif < 0L | Tarif > 5L){
    Tarif <- 1L
    warning("no inventory with indicator >5 (and < 0)! set inv=1 ")
  }

  if(identical(Tarif, 5L) & !requireNamespace("rBDAT", quietly = TRUE)){
    message("R-package 'rBDAT' not installed; please install.packages('rBDAT'); \n Tarif set to 0")
    Tarif <- 0L
  }

  if (identical(Tarif, 0L)){
    ## calculate the default taper from (quotient d03/d005) of taper curves
    d03 <- tprDiameter(tprTrees(spp=spp, Dm=Dbh, Hm=rep(1.3, length(Dbh)), Ht=Ht),
                       Hx = 0.3*Ht, cp = FALSE)
    d005 <- tprDiameter(tprTrees(spp=spp, Dm=Dbh, Hm=rep(1.3, length(Dbh)), Ht=Ht),
                        Hx = 0.05*Ht, cp = FALSE)
    EQ03 <- d03 / d005

  } else if(identical(Tarif, 5L)){
    if(requireNamespace("rBDAT", quietly = TRUE)){
      ## calculate the default taper from (quotient d03/d005) of taper curves
      EQ03 <- rBDAT::getForm(list(spp=spp, D1=Dbh, H1=1.3, H=Ht), inv=0)
    }

  } else {
    ## evaluate the taper form according to a specified inventory
    BDATSKNr <- BaMap(spp, type=9)

    a11 <- EQP[Tarif, 4, 1, BDATSKNr]
    a12 <- EQP[Tarif, 5, 1, BDATSKNr]
    a13 <- EQP[Tarif, 6, 1, BDATSKNr]

    h11 <- EQP[Tarif, 2, 1, BDATSKNr]
    h12 <- EQP[Tarif, 3, 1, BDATSKNr]
    h13 <- (h12 + h11) * 0.5

    D1 <- EQP[Tarif, 1, 1, BDATSKNr]

    a21 <- EQP[Tarif, 4, 2, BDATSKNr]
    a22 <- EQP[Tarif, 5, 2, BDATSKNr]
    a23 <- EQP[Tarif, 6, 2, BDATSKNr]

    h21 <- EQP[Tarif, 2, 2, BDATSKNr]
    h22 <- EQP[Tarif, 3, 2, BDATSKNr]
    h23 <- (h22 + h21) * 0.5

    D2 <- EQP[Tarif, 1, 2, BDATSKNr]

    Phi <- EQP[Tarif, 7, 1, BDATSKNr]

    ############################################################################
    # Z(H|D(i)) = MW [Q0.3| H | a(H(i,j)|D(i))]; j=1,2,3; i=1,2                #
    # Ratkowsky, D.A. (1990) (4.3.9), S97                                      #
    ############################################################################

    Q1 <- 2 * (Ht - h11) / (h12 - h11)
    Z1 <- a11 + (a12 - a11) *
      (1 - ((a12 - a13)/(a13 - a11))** Q1) /
      (1 - ((a12 - a13) / (a13 - a11)) ** 2)

    Q2 <- 2 * (Ht - h21) / (h22 - h21)
    Z2 <- a21 + (a22 - a21) *
      (1 - ((a22 - a23)/(a23 - a21)) ** Q2) /
      (1 - ((a22 - a23) / (a23 - a21)) ** 2)

    ############################################################################
    # EQ0.3(D,H) =  E [Q0.3| D, Z(H|D(i)); i=1,2]                              #
    # Ratkowsky, D.A. (1990) (4.3.23), S104                                    #
    ############################################################################

    EQ03 <- Z1 * Z2 * (D2 ** Phi - D1 ** Phi) /
      (Z2 * (D2 ** Phi - Dbh ** Phi) +
      Z1 * (Dbh ** Phi - D1 ** Phi))

    EQ03 <- ifelse(EQ03 < EQ03uG, EQ03uG, EQ03) # respect lower threshold
    EQ03 <- ifelse(EQ03 > EQ03oG, EQ03oG, EQ03) # respect upper threshold
  }
  return(EQ03)
}



getEQP <- function(){
  ## define size ---------------------------
  EQP <- array(dim=c(4, 7, 2, 9))
  # EQP[inv, param, set1/2, species]
  ## BWI1 ----------------------------------
  EQP[1, , , ] <- array(
    c(20, 10, 50, 0.650, 0.875, 0.850, 0.250,
      55, 10, 50, 0.525, 0.860, 0.775, 0.000,
      20, 10, 50, 0.750, 0.950, 0.860, 0.500,
      70, 10, 50, 0.670, 0.875, 0.780, 0.000,
      20, 10, 50, 0.650, 0.925, 0.860, 0.300,
      70, 10, 50, 0.500, 0.825, 0.700, 0.000,
      20, 10, 50, 0.700, 0.800, 0.770, 0.250,
      50, 10, 50, 0.700, 0.790, 0.770, 0.000,
      20, 10, 50, 0.725, 0.950, 0.875, 0.500,
      60, 10, 50, 0.630, 0.875, 0.750, 0.000,
      20, 10, 50, 0.700, 0.900, 0.830, 0.750,
      60, 10, 50, 0.650, 0.870, 0.820, 0.000,
      20, 10, 50, 0.700, 0.850, 0.840, 0.750,
      60, 10, 50, 0.675, 0.840, 0.825, 0.000,
      20, 10, 50, 0.775, 0.850, 0.810, 1.000,
      60, 10, 50, 0.725, 0.800, 0.760, 0.000,
      20, 10, 50, 0.700, 0.900, 0.830, 0.750,
      60, 10, 50, 0.650, 0.870, 0.820, 0.000),
    dim = c(7, 2, 9))
  # EQP[Inventur, Parameter, Set, Baumart]
  EQP[1, , 1, 1] # BWI1, alle parmeter, set1, Fichte
  EQP[1, , 2, 1] # BWI1, alle parmeter, set2, Fichte
  EQP[1, , 1, 9] # BWI1, alle parmeter, set1, Birke

  ## Neue Bundesländer ------------------------
  EQP[2, , ,] <- array(
    c(10, 10, 50, 0.700, 0.900, 0.875, 0.250,
      70, 10, 50, 0.550, 0.810, 0.760, 0.000,
      20, 10, 50, 0.750, 0.950, 0.860, 0.500,
      70, 10, 50, 0.670, 0.875, 0.780, 0.000,
      20, 10, 50, 0.650, 0.925, 0.860, 0.300,
      70, 10, 50, 0.500, 0.825, 0.700, 0.000,
      10, 10, 50, 0.710, 0.810, 0.780, 0.250,
      70, 10, 50, 0.700, 0.800, 0.760, 0.000,
      20, 10, 50, 0.725, 0.950, 0.875, 0.500,
      60, 10, 50, 0.630, 0.875, 0.750, 0.000,
      10, 10, 50, 0.625, 0.820, 0.800, 0.250,
      70, 10, 50, 0.575, 0.810, 0.780, 0.000,
      10, 10, 50, 0.650, 0.850, 0.825, 0.250,
      70, 10, 50, 0.575, 0.825, 0.800, 0.000,
      20, 10, 50, 0.775, 0.850, 0.810, 1.000,
      60, 10, 50, 0.725, 0.800, 0.760, 0.000,
      10, 10, 50, 0.700, 0.850, 0.825, 0.250,
      70, 10, 50, 0.575, 0.725, 0.675, 0.000),
    dim = c(7, 2, 9))
  # EQP[2, , 1, 1]

  ## Inventurstudie 2008 / IS08 ---------------------
  EQP[3, , ,] <- array(
    c(20, 10, 50, 0.494, 0.939, 0.857, 0.250,
      70, 10, 50, 0.494, 0.796, 0.729, 0.000,
      20, 10, 50, 0.665, 0.870, 0.838, 0.500,
      80, 10, 50, 0.593, 0.827, 0.745, 0.000,
      20, 10, 50, 0.519, 0.782, 0.768, 0.500,
      80, 10, 50, 0.467, 0.775, 0.755, 0.000,
      10, 10, 50, 0.674, 0.784, 0.770, 0.250,
      70, 10, 50, 0.650, 0.780, 0.761, 0.000,
      10, 10, 50, 0.689, 0.833, 0.813, 0.500,
      60, 10, 50, 0.417, 0.813, 0.753, 0.000,
      10, 10, 50, 0.467, 0.830, 0.806, 0.250,
      70, 10, 50, 0.400, 0.818, 0.781, 0.000,
      10, 10, 50, 0.586, 0.861, 0.827, 0.250,
      70, 10, 50, 0.482, 0.827, 0.793, 0.000,
      20, 10, 50, 0.775, 0.850, 0.810, 1.000,
      60, 10, 50, 0.725, 0.800, 0.760, 0.000,
      10, 10, 50, 0.454, 0.912, 0.812, 0.250,
      70, 10, 50, 0.502, 0.737, 0.680, 0.000),
    dim = c(7, 2, 9))
  # EQP[3,,1,7]

  ## BWI 3, 2012 / NFI 3 2012 _--------------
  EQP[4, , ,] <- array(
    c(20, 10, 50, 0.6420103, 0.888355, 0.8613915, 0.25,
      70, 10, 50, 0.4607479, 0.8679195, 0.7564586, 0.25,
      20, 10, 50, 0.7549, 0.9659322, 0.8680001, 0.5,
      80, 10, 50, 0.6573403, 0.8762392, 0.7838848, 0.5,
      20, 10, 50, 0.6397411, 0.8182462, 0.7964373, 0.5,
      80, 10, 50, 0.4664799, 0.7925813, 0.7392171, 0.5,
      10, 10, 50, 0.7204486, 0.967768, 0.7891799, 1,
      70, 10, 50, 0.7319144, 0.8527102, 0.7562992, 1,
      10, 10, 50, 0.788315, 0.8702361, 0.8701137, 0.5,
      60, 10, 50, 0.4612521, 0.8376158, 0.7608057, 0.5,
      10, 10, 50, 0.6985465, 0.8321258, 0.7983484, 0.25,
      70, 10, 50, 0.6078472, 0.8996977, 0.8125365, 0.25,
      10, 10, 50, 0.6925529, 0.8351468, 0.8111729, 0.25,
      70, 10, 50, 0.6565229, 0.8682255, 0.8281585, 0.25,
      20, 10, 50, 0.7684993, 0.9120466, 0.7992741, 1,
      60, 10, 50, 0.7262412, 0.8523916, 0.7576412, 1,
      10, 10, 50, 0.6772199, 0.962364, 0.8331439, 0.25,
      70, 10, 50, 0.5566832, 0.9492619, 0.7215377, 0.25),
    dim = c(7, 2, 9))
  # EQP[4,,1,7]
  return(EQP)
}

