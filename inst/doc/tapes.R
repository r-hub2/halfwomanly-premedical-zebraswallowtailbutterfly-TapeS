## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- eval=FALSE--------------------------------------------------------------
#  remotes::install_gitlab(„vochr/tapes“, build_vignettes = TRUE)

## ---- eval=FALSE--------------------------------------------------------------
#  help("tprTrees-class")

## ---- eval=FALSE--------------------------------------------------------------
#  help("parSort-class")

## -----------------------------------------------------------------------------
require(TapeS)
tprSpeciesCode(inSp = NULL, outSp = NULL)

## -----------------------------------------------------------------------------
BaMap(Ba = NULL, type = NULL)

## -----------------------------------------------------------------------------
require(TapeS)
obj <- tprTrees(spp=c(1, 3),
                Hm=list(c(1.3, 5, 7), c(1.3)),
                Dm=list(c(27, 23.5, 22.4), c(27)),
                Ht=c(27, 30))
Hx <- c(1.3, 5, 7)

## -----------------------------------------------------------------------------
tprDiameter(obj, Hx = Hx) # R0=TRUE, taper curve through measurements
tprDiameter(obj, Hx = Hx, bark = FALSE)
tprDiameter(obj, Hx = Hx, interval = "prediction")

## -----------------------------------------------------------------------------
tprDiameter(obj, Hx = 0.3*Ht(obj), interval = "prediction", cp=FALSE)

## -----------------------------------------------------------------------------
tprHeight(obj, Dx = c(10, 9, 8, 7))
tprHeight(obj, Dx = c(10, 9, 8, 7), bark = FALSE)

## -----------------------------------------------------------------------------
tprBark(obj, Hx = c(1, 2, 3))

## -----------------------------------------------------------------------------
tprVolume(obj) # default is Vfm
tprVolume(obj, AB = list(A=0, B=7), iAB=c("h", "dob"), bark=TRUE) # same
Vfm(obj) # wrapper
VolR(obj) # wrapper
Efm(obj, stH = 0.01) # default
VolE(obj)
VolFAO(obj)
Vfm_phys(obj) # takes a while
Efm_phys(obj)
tprVolume(obj, AB = list(A=0.01*Ht(obj), B=7, sl=0.01), iAB = c("H", "Dob"), bark=FALSE)

## -----------------------------------------------------------------------------
tprAssortment(obj) ## default assortment parameters
pars <- parSort(stH=0.2, Lxh=c(1, 1.5), fixN=2, fixL=4)
tprAssortment(obj, pars = pars)

## -----------------------------------------------------------------------------
tprBiomass(obj) # bwi-biomass

## -----------------------------------------------------------------------------
tprBiomass(obj, component = c("sw", "sb", "ndl")) 
tprBiomass(obj, component = c("all"))

## ---- eval=TRUE---------------------------------------------------------------
setTapeSoptions(Rfn = list(fn="zero"))
tprDiameter(obj, Hx=1.3)
setTapeSoptions(Rfn = list(fn="sig2"))
tprDiameter(obj, Hx=1.3)

## -----------------------------------------------------------------------------
plot(tprTrees(spp=3, Dm=7.9, Hm=1.3, Ht=12, ), mono=FALSE)

## -----------------------------------------------------------------------------
plot(tprTrees(spp=3, Dm=7.9, Hm=1.3, Ht=12, ), mono=TRUE)

