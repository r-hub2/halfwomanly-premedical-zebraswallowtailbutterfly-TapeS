#' @title total aboveground and component biomass
#' @description calculate total above ground and optionally component biomass
#' for given trees
#' @param obj object of class 'tprTrees'
#' @param component component for which biomass should be returned. If NULL,
#' total aboveground biomass is returned, if 'all', all components are returned.
#' See details.
#' @param useNFI if \code{TRUE}, agb is estimated by the NFI-functions and
#' component estimates are scaled so that their sum (i.e. agb) equals the
#' estimate of the NFI functions. If \code{FALSE}, the NSUR functions are used
#' for agb and component estimates.
#' @param interval character to indicate whether and which type of interval is
#' required; one of \code{none}, \code{confidence} or \code{prediction}.
#' @param mono logical, defaults to true. If calibrated taper curve is
#' non-monotonic at stem base, a support diameter is added.
#' @param Rfn Rfn setting for residuals error matrix, defaults to
#' \code{list(fn="sig2")}, see \code{\link[TapeR]{resVar}}.
#' @details The available components are agb (= total aboveground biomass),
#' stw (=stump wood), stb (=stump bark), sw (=solid wood with diameter above
#'  7cm over bark), sb (=bark of component sw), fwb (=fine wood incl. bark)
#'  and ndl (=needles), if applicable. The needles-component is set to zero for
#'  deciduous tree species, no mass for leaves is available. One can request
#'  'all' components to receive all components.
#' @return  a vector in case agb or only one component is requested, otherwise
#' a matrix with one row per tree
#' @references Kändler, G. and B. Bösch (2012). Methodenentwicklung für die
#' 3. Bundeswaldinventur: Modul 3 Überprüfung und Neukonzeption einer
#' Biomassefunktion - Abschlussbericht. Im Auftrag des Bundesministeriums für
#' Ernährung, Landwirtschaft und Verbraucherschutz in Zusammenarbeit mit dem
#' Institut für Waldökologie und Waldinventur des Johann Heinrich von
#' Thünen-Instituts, FVA-BW: 71.
#'
#' Kaendler (2021): Biometrische Modelle für die Ermittlung des Holzvorrats,
#' seiner Sortimentsstruktur und der oberirdischen Biomasse im Rahmen der
#' Bundeswaldinventur. Allg. Forst- u. J.-Ztg., 191. Jg., 5/6 83
#'
#' Vonderach, C., G. Kändler and C. Dormann (2018): Consistent set of additive
#' biomass equations for eight tree species in Germany fitted by nonlinear
#' seemingly unrelated regression. Annals of Forest Science (2018) 75:49
#' doi: 10.1007/s13595-018-0728-4

#' @export
#' @examples
#' obj <- tprTrees(spp=c(1, 15),
#'                 Dm=list(c(30, 28), c(30, 28)),
#'                 Hm=list(c(1, 3), c(1, 3)),
#'                 Ht = rep(30, 2))
#' (tmp <- tprBiomass(obj, component="all"))
#'
#' tprBiomass(obj, component=NULL) # aboveground biomass
#' component <- c("agb", "sw", "sb", "ndl")
#' tprBiomass(obj, component=component)
#' component <- c("sw", "sb", "ndl")
#' tprBiomass(obj, component="all")
#' # use NSUR-functions from Vonderach et al. 2018
#' # obs: currently sth=1% of tree height
#' # and kl=70% of tree height
#' tprBiomass(obj, component="all",  useNFI = FALSE)
#'
#' ## getting confidence and prediction intervals
#' useNFI <- FALSE
#' interval <- "confidence"
#' component <- c("sw", "agb")
#' tprBiomass(obj, component, useNFI, interval)
#' tprBiomass(obj, component, useNFI, interval="none")
#' tprBiomass(obj, component, useNFI=TRUE, interval)
#' tprBiomass(obj, component, useNFI=TRUE, interval="none")
#'
#' obj <- tprTrees(spp=15, Dm=30, Hm=1.3, Ht=27)
#' tprBiomass(obj, component="all", interval="confidence")
#' tprBiomass(obj, component="ndl", interval="confidence")
#'
#' obj <- tprTrees(spp=c(1, 15), Dm=c(30, 30), Hm=c(1.3, 1.3), Ht=c(27, 27))
#' tprBiomass(obj, component="all", interval="confidence")
#'
#' obj <- tprTrees(spp=c(1, 15), Dm=c(30, 30), Hm=c(1.3, 1.3), Ht=c(27, 27))
#' tprBiomass(obj, component=c("sw", "ndl"), interval="confidence")
#'
#' obj <- tprTrees(spp=c(1, 15), Dm=c(30, 30), Hm=c(1.3, 1.3), Ht=c(27, 27))
#' tprBiomass(obj, component=c("ndl", "agb"), interval="confidence")
#'
#' obj <- tprTrees(spp=c(1, 15), Dm=c(30, 30), Hm=c(1.3, 1.3), Ht=c(27, 27))
#' tprBiomass(obj, component=c("ndl"), interval="confidence")


setGeneric("tprBiomass",
           function(obj, component=NULL, useNFI=TRUE, interval="none", mono=TRUE, Rfn=NULL)
             standardGeneric("tprBiomass"))

#' @describeIn tprBiomass method for class 'tprTrees'
setMethod("tprBiomass", signature = "tprTrees",
          function(obj, component=NULL, useNFI=TRUE, interval="none", mono=TRUE, Rfn=NULL){

            if(!is.null(Rfn)){
              old <- options()$TapeS_Rfn
              options(TapeS_Rfn = Rfn)
              on.exit(options(TapeS_Rfn = old))
            }
            if(!identical(mono, oldmono <- options()$TapeS_mono)){
              options("TapeS_mono" = mono)
              on.exit(options("TapeS_mono" <- oldmono))
            }
            if(!(interval[1] %in% c("none", "confidence", "prediction"))){
              interval <- "none"
            } else {
              interval <- interval[1]
            }

            ## check component names on validity
            # check_Comp <- TapeS:::check_Comp
            comp <- check_Comp(component)

            d13 <- D13(obj)
            d03 <- D03(obj)


            if(isTRUE(useNFI) & identical(comp, "agb") & identical(interval, "none")){

              ## getting the agb biomass from Riedel & Kaendler 2017
              bwiBm <- biomass(BaMap(obj@spp, 6), d13, d03, obj@Ht)
              return(bwiBm)

            } else {

              ## either several components OR raw NSUR estimates are required
              ## hence, we need to evaluate the NSUR functions
              ## get the component shares from NSUR-Cpp-Functions
              ## Vonderach et al. 2018 (but re-fitted using D03 and KL)
              sppBmC <- BaMap(obj@spp, 7) # ba code mapping to comp-functions
              nsurBm <- as.data.frame( nsur(sppBmC, dbh = d13, ht = obj@Ht,
                                            sth = 0.01*obj@Ht, d03 = d03,
                                            kl = 0.7 * obj@Ht) )
              nsurBm$agb <- rowSums(nsurBm[, -which(colnames(nsurBm)=="id")])

              if(interval %in% c("confidence", "prediction")){ # implementing variance estimate
                # NSURvar <- TapeS:::NSURvar
                # comp <- "all"
                # comp <- TapeS:::check_Comp(comp)
                vcomp <- NSURvar(data = data.frame(spp = sppBmC, dbh = d13, ht = obj@Ht,
                                                   sth = 0.01*obj@Ht, D03 = d03,
                                                   kl = 0.7 * obj@Ht),
                                 estBM = nsurBm[, comp, drop=FALSE], comp = NULL,
                                 interval = interval, level=0.95, adjVarPar=TRUE, as.list=TRUE)
              } else {
                vcomp <- NULL
              }

              if(isTRUE(useNFI)){ # NFI estimate required
                bmshare <- nsurBm[, comp, drop=FALSE] / nsurBm$agb
                bwiBm <- biomass(BaMap(obj@spp, 6), d13, d03, obj@Ht)
                res <- apply(bmshare[, comp, drop=FALSE],
                             MARGIN = 2,
                             function(m){ bwiBm * m }, simplify = FALSE)
                res <- as.data.frame(res)

                if(!is.null(vcomp)){
                  # vcomp holds the upper and lower bound and MSE estimate
                  for(i in seq(along=vcomp)){ # iterate over components
                    # i <- 1
                    nm <- names(vcomp[i])
                    lwr <- paste0(nm, "_lwr")
                    upr <- paste0(nm, "_upr")
                    ECBM <- paste0(nm, "_ECBM")
                    MSE <- paste0(nm, "_MSE")

                    # proportionally re-scale MSE and interval bounds
                    vc <- vcomp[[i]]
                    vc$Q <- vc[ , upr] / vc[, ECBM] - 1
                    vc$qt <- (vc[, upr] - vc[, ECBM]) / sqrt(vc[, MSE])
                    vc[, ECBM] <- res[, nm]
                    vc[, lwr] <- vc[, ECBM] * (1 - vc$Q)
                    vc[, upr] <- vc[, ECBM] * (1 + vc$Q)
                    vc[, MSE] <- ((vc[, upr] - vc[, ECBM]) / vc$qt)^2
                    vc$Q <- vc$qt <- NULL
                    vc[which(vcomp[[i]]==0, arr.ind = TRUE)] <- 0
                    vcomp[[i]] <- vc
                  }
                  names(vcomp) <- NULL
                  res <- do.call(cbind, vcomp)
                  colnames(res) <- gsub("_ECBM$", "", colnames(res))
                  return(res)

                } else {
                  # vcomp==NULL => only point estimate required
                  return(res)
                }


              } else { # NSUR estimate required

                if(!is.null(vcomp)){
                  # vcomp holds the upper and lower bound and MSE estimate
                  names(vcomp) <- NULL
                  res <- do.call(cbind, vcomp)
                  colnames(res) <- gsub("_ECBM$", "", colnames(res))
                  return(res)
                } else {
                  # vcomp==NULL => only point estimate required
                  res <- nsurBm[, comp]
                  return(res)
                }

              }
            }
          })

