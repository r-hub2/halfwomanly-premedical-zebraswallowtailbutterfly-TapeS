#' function for mapping the 36 tree species to several internal functions
#' @param Ba BDAT tree number
#' @param type a number referring to the type to be mapped
#' @details c.f. BDAT source code, line 7622, data block Ban(36, 7)
#' type 1: Schaftform // taper form
#' type 2: Rinde // bark
#' type 3: Durchschnittliche Aufarbeitungsgrenze (nach EST) //average cutting diameter
#' type 4: Höhe unverwertbares Derbholz // percentage non-merchantable coarse wood
#' type 5: durchschnittlicher Astdurchmesser in der Krone // average branch diameter inside crown
#' type 6: BWI-Biomasse-Funktionen // NFI-biomass functions according to Riedel & Kändler (2017)
#' type 7: kompartimentweise Biomassefunktionen // component biomass functions according to Vonderach et al (2018)
#' type 8: Zuordnung zu volfao // Mapping to volume according to FAO (FIX: mapping still temporary)
#' Not included: volume tables according to Grundner and Schwappach as well as
#' volume tables according to Krenn for small trees below 10cm dbh
#' @return value(s), either a scalar, vector or matrix, with respect to tree
#' species mapping to functions
#' @export
#' @examples
#' BaMap(1,1) # which taper form for Norway spruce
#' BaMap(15,1) # which taper form for European Beech
#' BaMap(15,2) # which bark equation for European Beech
#' BaMap(,1) # return all taper form mappings
#' BaMap(1,) # return all mappings for Norway spruce
#' BaMap() # return all mappings
#' BaMap(, 6) # biomass mapping
#' BaMap(, 7) # component biomass functions
#' BaMap(, 8) # mapping for Vol_FAO
BaMap <- function(Ba=NULL, type=NULL){

  if(!is.null(type)){
    if(any(type > 9) | any(type < 1)){
      stop("'type' must be integer inside interval [1, 9]!")
    }
  } else {
    type <- 1:9
  }
  if(!is.null(Ba)){
    if(any(Ba > 39) | any(Ba < 1)){
      stop("'Ba' must be integer inside interval [1,39]!")
    }
  } else {
    Ba <- 1:39
  }

  mapBa <- array(NA, dim = c(39, 9),
                 dimnames = list(1:39, c("Schaftform", "Rinde", "Az", "uvDh",
                                         "mAS", "Biomasse", "BmComp", "VolFao",
                                         "Formigkeit")))
  ## Zuordnung der Schaftform
  mapBa[, 1] <- c(1, 1, 2, 2, 4, 4, 4, 3, 5, 5, 5, 1, 1, 1, 6, 6, 7, 8,
                  8, 8, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 6, 6, 6,
                  9, 10, 11)
  ## Zuordnung der Rindengleichungen
  mapBa[, 2] <- c(1, 1, 2, 2, 4, 8, 9, 3, 6, 6, 7, 1, 1, 2,10,18,11,14,
                  25,25,17,15,15,15,21,21,16,28,26,20,19,23,11,12,12,26,
                  10, 11, 14)
  ## Zuordnung der durchschnittlichen Aufarbeitungsgrenze
  mapBa[, 3] <- c(1, 1, 3, 3, 4, 4, 4, 1, 5, 5, 5, 1, 1, 1, 6, 6, 7, 7,
                  6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 6, 6, 6, 6, 7,
                  6, 7, 7)
  ## Zuordnung der Höhe des unverwertbaren Derbholzes in der Krone (Laubholz)
  mapBa[, 4] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 2,
                  2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2,
                  1, 2, 2)
  ## Zuordnung der durchschnittlichen Astdurchmesser in der Krone (Laubholz)
  mapBa[, 5] <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 2,
                  2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1,
                  1, 2, 2)
  ## Zuordnung zur Biomasse-Baumart, vgl. Koeff.f, Zeile 133f
  mapBa[, 6] <- c(1,1,2,2,4,4,4,3,5,5,5,1,1,1,6,9,7,7,17,17,8,10,10,11,
                  10,15,12,16,6,14,13,6,6,18,6,6,
                  6, 7, 7)
  ## Zuordnung zur Biomassekompartiment-Baumart
  ## 1=fi, 2=ta, 3=kie, 4=dgl, 5=bu, 6=ei, 7=bah, 8=es
  mapBa[, 7] <- c(1,1,2,2,3,3,3,4,4,4,4,1,1,1,5,6,6,6,5,5,8,7,
                  7,7,7,5,5,5,5,5,5,5,5,5,5,5,
                  5, 6, 6)
  ## Zuordnung zur Tabelle volfao, FIX: this is temporary
  mapBa[, 8] <- c(1,1,2,2,4,4,4,3,5,5,5,1,1,1,6,7,7,7,6,6,7,7,
                  7,7,7,6,6,6,6,6,6,6,6,6,6,6,
                  6, 7, 7)
  ## Zuordnung zur Formigkeit verschiedener Zeitpunkte
  ## 1=Fi, 2=Ta, 3=Dgl, 4=Kie, 5=Lae, 6=Bu, 7=Eiche, 8=Roteiche, 9=Birke
  ## vgl. BDAT Quellcode Formtarife.f
  mapBa[, 9] <- c(1, 1, 2, 2, 4, 4, 4, 3, 5, 5, 5, 1, 1, 1,
                  6, 6, 7, 8, 8, 8, 6, 6, 6, 6, 6, 9, 6, 6, 6, 6, 6, 6, 7, 6, 6, 6,
                  6, 7, 8)
  return(mapBa[Ba, type, drop=TRUE])
}

