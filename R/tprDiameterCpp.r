#' @title Function to extract diameters from Taper curve using Rcpp
#' @description This function uses Rcpp and C-code to implement the diameter
#' estimation of package TapeR to allow for faster estimation if no interval
#' information is required.
#' @param obj object of class 'tprTrees'
#' @param Hx vector of heights for which diameter are required
#' @param bark should diameter over or under bark be returned?
#' @param cp cartesian product, i.e. apply all \code{Hx} to all trees, defaults
#' to TRUE
#' @param mono logical to decide whether a supporting diameter should be added
#' in case the taper curve is regarded as non-monotonic. Defaults to TRUE.
#' @param Rfn setting for residuals error matrix, defaults to \code{"sig2"}, see
#' details.
#' @details Function evaluates taper curves at required height \code{Hx}. By
#' default (\code{cp==TRUE}), the taper curve is evaluated at \code{Hx} for each
#' tree. If \code{cp==FALSE}, each tree is evaluated at exactly one Hx (recycled
#' if necessary). This feature is intended for situations where diameter in
#' relative heights are required. Then, the recycling of one height Hx (e.g.
#' 1.3m) is not possible, since relative heights depend on absolute tree height,
#' which might be different for each tree. Hence a call like
#' \code{tprDiameter(obj, Hx=0.3*Ht(obj), cp=FALSE)} is necessary.
#' @return a vector, in case only one diameter (i.e. Hx) is required per tree
#' (\code{cp=FALSE}) or a matrix of size
#' \code{length(trees)} x \code{length(Hx)} (\code{cp=TRUE}).
#' @seealso \code{\link{tprDiameter}} if confidence or prediction intervals
#' are required.
#' @import TapeR
#' @import Rcpp
#' @export
#' @useDynLib TapeS, .registration = TRUE
#' @examples
#' obj <- tprTrees(spp=c(1 , 3),
#'                 Hm=list(c(1.3, 5), c(1.3, 5)),
#'                 Dm=list(c(27, 25), c(27, 25)),
#'                 Ht=c(27, 27))
#' Hx <- seq(0, 1, 0.1)
#' tprDiameterCpp(obj, Hx = Hx)
#' tprDiameterCpp(obj, Hx = Hx, bark=FALSE)
#' tprDiameterCpp(obj, Hx = c(1, 2), bark=FALSE, cp=FALSE)
#'
#' \donttest{
#' require(rbenchmark)
#' benchmark(tprDiameter(obj, Hx, bark = TRUE),
#'           tprDiameterCpp(obj, Hx, bark = TRUE),
#'           replications = 10000)[,1:4]
#' }

setGeneric("tprDiameterCpp",
           function(obj, Hx, bark=TRUE, cp=TRUE, mono=TRUE, Rfn=NULL)
             standardGeneric("tprDiameterCpp"))

#' @describeIn tprDiameterCpp method for class 'tprTrees'
setMethod("tprDiameterCpp", signature = "tprTrees",
          function(obj, Hx, bark=TRUE, cp=TRUE, mono=TRUE, Rfn=NULL){
            SKspp <- BaMap(obj@spp, 1)
            if(is.null(Rfn)) Rfn <- getOption("TapeS_Rfn")
            if(isTRUE(cp)){
              ## cartesian product
              ## apply all Hx to all trees subsequently via Hx[TRUE] => all Hx
              use <- rep(TRUE, length(obj))
            } else {
              ## no cartesian product between obj and Hx
              ## for each tree, use exactly one given Hx
              if(!identical(length(obj), length(Hx))){
                stop("if cp='FALSE', length of 'obj' and 'Hx' must be identical.")
              } else {
                use <- seq(along=obj)
              }
            }
            ncol <- ifelse(all(use==TRUE), length(Hx), 1)
            res <- matrix(NA, nrow = length(obj@spp), ncol = ncol)
            if(isTRUE(bark)){
              for(i in seq(along=obj@spp)){
                if(isFALSE(obj@monotone[i]) & isTRUE(mono)){
                  suppD <- attr(obj@monotone, "D005")[i] * SKPar[[ SKspp[i] ]]$q001
                  suppH <- 0.01*obj@Ht[i]
                } else {
                  suppD <- suppH <- numeric(0)
                }
                rv <- TapeR::resVar(c(suppH, obj@Hm[[i]])/obj@Ht[i], fn=Rfn$fn,
                                     sig2 = SKPar[[ SKspp[i] ]]$sig2_eps,
                                     par = Rfn$par)
                ## calling c++ function lmeSKEBLUP
                res[i, ] <- lmeSKEBLUP(xm=c(suppH, obj@Hm[[i]])/obj@Ht[i],
                                       ym=c(suppD, obj@Dm[[i]]),
                                       xp=Hx[ use[i] ]/obj@Ht[i],
                                       par=SKPar[[ SKspp[i] ]],
                                       RV=rv)$yp
              }
            } else {
              for(i in seq(along=obj@spp)){
                if(isFALSE(obj@monotone[i]) & isTRUE(mono)){
                  suppD <- attr(obj@monotone, "D005")[i] * SKPar[[ SKspp[i] ]]$q001
                  suppH <- 0.01*obj@Ht[i]
                } else {
                  suppD <- suppH <- numeric(0)
                }
                rv <- TapeR::resVar(c(suppH, obj@Hm[[i]])/obj@Ht[i], fn=Rfn$fn,
                                     sig2 = SKPar[[ SKspp[i] ]]$sig2_eps,
                                     par = Rfn$par)
                ## calling c++ function lmeSKEBLUP
                DHx <- lmeSKEBLUP(xm=c(suppH, obj@Hm[[i]])/obj@Ht[i],
                                  ym=c(suppD, obj@Dm[[i]]),
                                  xp=Hx[ use[i] ]/obj@Ht[i],
                                  par=SKPar[[ SKspp[i] ]],
                                  RV=rv)$yp
                res[i, ] <- DHx - bark(rep(obj@spp[i], length(DHx)), DHx,
                                       relH = Hx[ use[i] ]/obj@Ht[i])
              }
            }
            return(res[,, drop=TRUE])
          })

