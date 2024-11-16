#' @title Estimate height for given diameter w/ or w/o bark
#' @description Function to extract the height of given diameter w/ or w/o bark
#' from taper curve
#' @param obj object of class 'tprTrees'
#' @param Dx diameter for which height is required
#' @param bark should given diameter be considered over or under bark?
#' @param cp cartesian product, i.e. apply all \code{Hx} to all trees, defaults
#' to TRUE
#' @param mono logical, defaults to true. If calibrated taper curve is
#' non-monotonic at stem base, a support diameter is added.
#' @param Rfn Rfn setting for residuals error matrix, defaults to
#' \code{list(fn="sig2")}, see \code{\link[TapeR]{resVar}}.
#' @return estimated height of given diameter
#' @seealso \code{\link{tprDiameter}}, \code{\link{tprDiameterCpp}}
#' @import TapeR
#' @export
#' @examples
#' obj <- tprTrees(spp=c(1, 3, 8, 15),
#'                 Dm=list(c(30, 28), c(30, 28), c(30, 28), c(30, 28)),
#'                 Hm=list(c(1.3, 5), c(1.3, 5), c(1.3, 5), c(1.3, 5)),
#'                 Ht = rep(30, 4))
#' tprHeight(obj, Dx = c(30, 7), bark=TRUE)
#' tprHeight(obj, Dx = c(30, 7), bark=FALSE)
#'
#' ## no cartesion product between obj and Dx, i.e. cp=FALSE
#' ## Dx is recycled if necessary
#' tprHeight(obj, Dx = c(30, 7), bark=FALSE, cp=FALSE)

setGeneric("tprHeight",
           function(obj, Dx, bark=TRUE, cp=TRUE, mono=TRUE, Rfn=NULL)
             standardGeneric("tprHeight"))

#' @describeIn tprHeight method for class 'tprTrees'
setMethod("tprHeight", signature = "tprTrees",
          function(obj, Dx, bark=TRUE, cp=TRUE, mono=TRUE, Rfn=NULL){
            SKspp <- BaMap(obj@spp, 1)
            if(is.null(Rfn)) Rfn <- getOption("TapeS_Rfn")

            if(isTRUE(cp)){
              ## apply all Dx to all trees subsequently via Hx[TRUE] => all Hx
              use <- rep(TRUE, length(obj))
            } else {
              ## for each tree, use exactly one given Dx
              use <- seq(along=obj)
              # recycle Hx to be at least of obj-length
              Dx <- Dx * rep(1, max(length(obj)*length(Dx), length(Dx)))
            }
            ncol <- ifelse(all(use==TRUE), length(Dx), 1)
            res <- matrix(NA, nrow = length(obj@spp), ncol = ncol)
            for(i in seq(along=obj@spp)){
              if(isFALSE(obj@monotone[i]) & isTRUE(mono)){
                suppD <- attr(obj@monotone, "D005")[i] * SKPar[[ SKspp[i] ]]$q001
                suppH <- 0.01*obj@Ht[i]
              } else {
                suppD <- suppH <- numeric(0)
              }
              if(isFALSE(bark)){
                res[i, ] <- sapply(seq(along=Dx[ use[i] ]), function(a){
                  E_HDxoR_HmDm_Ht.f(DxoR = Dx[ use[i] ][a],
                                    Dm = c(suppD, obj@Dm[[i]]),
                                    Hm = c(suppH, obj@Hm[[i]]),
                                    mHt = obj@Ht[i],
                                    sHt = 0,
                                    par.lme = SKPar[[ SKspp[i] ]],
                                    Rfn = Rfn)
                })
              } else {
                res[i, ] <- sapply(seq(along=Dx[ use[i] ]), function(a){
                  TapeR::E_HDx_HmDm_HT.f(Dx = Dx[ use[i] ][a],
                                         Dm = c(suppD, obj@Dm[[i]]),
                                         Hm = c(suppH, obj@Hm[[i]]),
                                         mHt = obj@Ht[i],
                                         sHt = 0,
                                         par.lme = SKPar[[ SKspp[i] ]],
                                         Rfn = Rfn)
                })
              }
            }
            return(res)
          })
