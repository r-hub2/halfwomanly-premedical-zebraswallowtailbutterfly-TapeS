#' @title tree volume information
#' @description Wrapper to get specific type of volume from taper curve
#' @param obj a object of class 'tprTrees'
#' @details wrapper functions around \code{\link{tprVolume}}, which return
#' specific definitions of stem volume.
#' @return vector of volume estimates
#' @export
#' @examples
#' t <- tprTrees() # constructor of class 'tprTrees'
#' Vfm(t)
#' Efm(t)
#' Efm(t, stH=0.01) # stump height = 1\% of tree height
#' Efm(t, stH=10) # stump height=10cm
#' VolR(t)
#' VolE(t)
#' VolFAO(t)
#' Vfm_phys(t) # slower since much more evaluations of taper curve (every 1 cm)
#' Efm_phys(t, stH=0.01) # slower since much more evaluations of taper curve (every 1 cm)
Vfm <- function(obj){
  Vol <- tprVolume(obj, AB = list(A=0, B=7, sl=2), iAB = c("H", "Dob"), bark=TRUE)
  return(Vol)
}

#' @describeIn Vfm Efm, i.e. coarse wood excl. bark from Ht=stH*Ht to Dob=7cm
#' @param stH assumed or known relative or absolute stump height, from which
#' volume calculation should starts, defaults to 0.01
#' @details Function \code{Efm} uses parameter \code{stH} to define starting
#' point, i.e. stump height, of volume calculation. \code{stH} can be defined
#' relative to total tree height \code{(0 < stH <= 1)} or in absolute measure
#' (unit=cm) in case \code{stH > 1}
#'
#' @export
Efm <- function(obj, stH=0.01){
  stopifnot(is.numeric(stH), stH>0, length(stH) == 1 | length(stH) == length(obj@Ht))
  if(length(stH) != length(obj@Ht)){
    # then length(stH)==1 & length(obj@Ht)>1 and vector should be extended
    stH <- rep(stH, length(obj@Ht))
  }
  A <- ifelse(stH < 1, stH * obj@Ht, stH / 100) # stH in absolute meter [m]
  stopifnot(all(A < obj@Ht))
  Vol <- tprVolume(obj, AB = list(A=A, B=7, sl=2), iAB = c("H", "Dob"), bark=FALSE)
  return(Vol)
}

#' @describeIn Vfm VolR: Volume from H=0 to D=7cm over bark, measured as 2m sections
#' @export
VolR <- function(obj){
  Vol <- tprVolume(obj, AB = list(A=0, B=7, sl=2), iAB = c("H", "Dob"), bark=TRUE)
  return(Vol)
}

#' @describeIn Vfm VolE: sum of volume of default assortments according to RVR
#' @details \code{VolE} calculates as the sum of volume of default assortments
#'  (stem wood, top log, industrial wood, X-wood, non-usuable wood according to
#'  RVR. For dbh < 7cm a linear regression is applied.
#' @export
VolE <- function(obj){
  Vol <- tprAssortment(obj)
  Vol <- aggregate(vol ~ tree, data=Vol, FUN = "sum")
  dz <- Dbh(obj)*10
  vold7 <- 0.000008951 * dz^2 - 0.000850412 * dz + 0.018412163
  Vol$vol <- ifelse(dz < 70, vold7, Vol$vol)
  return(Vol)
}

#' @describeIn Vfm VolFAO: from stump to tree top incl. bark; if dbh < 7cm using
#' tabulated values
#' @details \code{VolFAO} calculates tree volume starting from stump up to tree
#' top (in contrast to german definition, which uses D=7cm over bark), and
#' includes bark component. Stump height is defined as 1\% of tree height.
#' Volume calculation is based on 2m-sections. For trees with dbh < 7cm,
#' tabulated values are used, see Riedel et al. (2017) for details (e.g. p.35,
#' table 5.6).
#' @references \href{https://www.bundeswaldinventur.de/fileadmin/SITE_MASTER/content/Downloads/BWI_Methodenband_web.pdf}{Riedel, T. and Hennig, P. and Kroiher, F. and Polley, H. and
#' Schwitzgebel, F. (2017): Die dritte Bundeswaldinventur (BWI 2012). Inventur-
#' und Auswertemethoden. 124 pages.}
#' @export
VolFAO <- function(obj){
  Vol <- tprVolume(obj, AB = list(A=0.01*obj@Ht, B=obj@Ht, sl=2), iAB = "H", bark=TRUE)
  dz <- Dbh(obj)*10
  dbhclass <- ifelse(dz < 1, 2, # dbh=0
                     ifelse(dz < 50, 3, ifelse(dz < 60, 4, 5)))
  sp <- BaMap(obj@spp, 8)
  vold7 <- sapply(1:length(sp), function(a) volfaodlt7[ sp[a], dbhclass[a] ])
  Vol <- ifelse(dz < 70, vold7, Vol)
  return(Vol)
}

#' @describeIn Vfm Vfm_phys physical volume of tree incl. bark from A=0
#' @details \code{Vfm_phys} is equal to \code{Vfm}, except that the taper curve
#' is numerically integrated, by use of section length of 0.01m. This is
#' relevant if biomass or nutrient export is to be calculate. Numerical
#' integration is quite slow.
#' @export
Vfm_phys <- function(obj){
  Vol <- tprVolume(obj, AB = list(A=0, B=7, sl=.01), iAB = c("H", "Dob"), bark=TRUE)
  return(Vol)
}

#' @describeIn Vfm Efm_phys physical volume of tree excl. bark from A=0.1*Ht
#' @details \code{Efm_phys} is equal to \code{Efm}, except that the taper curve
#' is numerically integrated, by use of section length of 0.01m. This is
#' relevant if biomass or nutrient export is to be calculate. Numerical
#' integration is quite slow.
#' @export
Efm_phys <- function(obj, stH=0.01){
  stopifnot(is.numeric(stH), stH>0, length(stH) == 1 | length(stH) == length(obj@Ht))
  if(length(stH) != length(obj@Ht)){
    # then length(stH)==1 & length(obj@Ht)>1 and vector should be extended
    stH <- rep(stH, length(obj@Ht))
  }
  A <- ifelse(stH < 1, stH * obj@Ht, stH / 100) # stH in absolute meter [m]
  stopifnot(all(A < obj@Ht))
  Vol <- tprVolume(obj, AB = list(A=A, B=7, sl=0.01), iAB = c("H", "Dob"), bark=FALSE)
  return(Vol)
}
