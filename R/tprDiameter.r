#' @title Functions to extract diameters from Taper curve
#' @description Function evaluates TapeR taper curve models for given trees
#' according to species, required height and optionally substracts double bark
#' thickness.
#' @param obj object of class 'tprTrees'
#' @param Hx vector of heights for which diameter w/ or w/0 bark are required
#' @param bark should diameter over or under bark be returned?
#' @param interval indicator about whether 'confidence' or 'prediction'
#' intervals are required (defaults to 'none'), optionally function returns the
#' mean squared error of the mean and predictions ('MSE').
#' @param cp cartesian product, i.e. apply all \code{Hx} to all trees, defaults
#' to TRUE
#' @param mono logical, defaults to true. If calibrated taper curve is
#' non-monotonic at stem base, a support diameter is added.
#' @param Rfn Rfn setting for residuals error matrix, defaults to
#' \code{list(fn="sig2")}, see \code{\link[TapeR]{resVar}}.
#' @details Function evaluates taper curves at required height \code{Hx}. By
#' default (\code{cp==TRUE}), the taper curve is evaluated at \code{Hx} for each
#' tree. If \code{cp==FALSE}, each tree is evaluated at exactly one Hx (recycled
#' if necessary). This feature is intended for situations where diameter in
#' relative heights are required. Then, the recycling of one height Hx (e.g.
#' 1.3m) is not possible, since relative heights depend on absolute tree height,
#' which might be different for each tree. Hence a call like
#' \code{tprDiameter(obj, Hx=0.3*Ht(obj), cp=FALSE)} is necessary.
#'
#' @return a matrix or data.frame depending on value of \code{interval}. If
#' 'none' (the default), a matrix of size [length(obj@Ht), length(Hx)] is
#' returned, otherwise a data.frame of size [length(obj@Ht) * length(Hx), 5].
#' The five columns hold a tree identifier, Hx, lower confidence/prediction
#' interval, the estimated diameter and the upper confidence/prediction
#' interval. In case 'interval=MSE' the returned columns contain a tree
#' identifier, Hx, the estimated diameter and mean squared error (MSE) of the
#' mean and of the prediction. Estimates and intervals include bark or not,
#' depending on \code{bark}.
#' @seealso \code{\link{tprDiameterCpp}} for a faster implementation if no
#' confidence or prediction information are required and \code{\link{tprBark}}
#' for the applied bark reduction.
#' @import TapeR
#' @export
#' @examples
#' ## prediction for new tree using implemented 'TapeR' taper curve model
#' obj <- tprTrees(spp=c(1, 3),
#'                 Hm=list(c(1.3, 5), c(1.3, 5)),
#'                 Dm=list(c(27, 25), c(27, 25)),
#'                 Ht=c(27, 27))
#' hx <- c(1.3, 5, 7)
#' ## by default, Hx applied on each tree, i.e. result is a 2x3 matrix
#' tprDiameter(obj, Hx = hx)
#'
#' ## if cp=FALSE, each tree only 'sees' one Hx, i.e. results is a vector
#' ## (obs: length of Hx must be identical to length of obj)
#' tprDiameter(obj, Hx = c(1.3, 5), cp=FALSE)
#' tprDiameter(obj, Hx = hx, bark = FALSE)
#' tprDiameter(obj, Hx = hx, interval = "confidence")
#' tprDiameter(obj, Hx = hx, bark = FALSE, interval = "prediction")
#' tprDiameter(obj, Hx = hx, interval = "MSE")
#' tprDiameter(obj, Hx = hx, bark=FALSE, interval = "MSE")
#'
#' ## here same behaviour, if cp=FALSE
#' tprDiameter(obj, Hx = c(1.3, 5), bark = FALSE,
#'             interval = "prediction", cp=FALSE)
#' ## using Cpp-implementation
#' ## faster, but no intervals available
#' tprDiameterCpp(obj, Hx = hx)
#' tprDiameterCpp(obj, Hx = c(1.3, 5), cp=FALSE)
#'
#' ## prediction for objects of class 'datBDAT':
#' if(require(rBDAT)){
#'   tree <- rBDAT::buildTree(list(spp=1, D1=20:30, H1=1.3, H2=50, H=20:30))
#'   tree <- bdat_as_tprtrees(tree)
#'   tprDiameter(tree, Hx = 1.3)
#' }

setGeneric("tprDiameter",
           function(obj, Hx, bark=TRUE, interval="none", cp=TRUE, mono=TRUE, Rfn=NULL)
             standardGeneric("tprDiameter"))

#' @describeIn tprDiameter method for class 'tprTrees'
setMethod("tprDiameter", signature = "tprTrees",
          function(obj, Hx, bark=TRUE, interval="none", cp=TRUE, mono=TRUE, Rfn=NULL){
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
            if(interval %in% c("confidence", "prediction")){
              ncol <- 5
              nrow <- 0
              colnms <- c("tree", "Hx", "lwr", "EDx", "upr")
            } else if(interval == "MSE"){
              ncol <- 5
              nrow <- 0
              colnms <- c("tree", "Hx", "EDx", "MSE_Mean", "MSE_Pred")
            } else if(interval == "none"){
              ncol <- ifelse(all(use==TRUE), length(Hx), 1)
              nrow <- length(obj)
              colnms <- NULL
            } else {
              stop("'interval' must be one of 'none', 'confidence', 'prediction' or 'MSE'")
            }
            res <- matrix(NA, ncol = ncol, nrow = nrow)
            if(isTRUE(bark)){
              for(i in seq(along=obj@spp)){
                if(isFALSE(obj@monotone[i]) & isTRUE(mono)){
                  suppD <- attr(obj@monotone, "D005")[i] * SKPar[[ SKspp[i] ]]$q001
                  suppH <- 0.01*obj@Ht[i]
                } else {
                  suppD <- suppH <- numeric(0)
                }
                tmp <- TapeR::E_DHx_HmDm_HT.f(Hx=Hx[ use[i] ],
                                              Dm = c(suppD, obj@Dm[[i]]),
                                              Hm = c(suppH, obj@Hm[[i]]),
                                              mHt = obj@Ht[i],
                                              sHt = obj@sHt[i],
                                              par.lme = SKPar[[ SKspp[i] ]],
                                              Rfn = Rfn)
                if(interval=="none"){
                  res[i, ] <- tmp$DHx
                } else if(interval=="confidence"){
                  res <- rbind(res, cbind(tree=i, Hx=Hx[ use[i] ], tmp$CI_Mean))
                } else if(interval=="prediction"){
                  res <- rbind(res, cbind(tree=i, Hx=Hx[ use[i] ], tmp$CI_Pred))
                } else if(interval=="MSE"){
                  res <- rbind(res, cbind(tree=i, Hx=Hx[ use[i] ], tmp$DHx, tmp$MSE_Mean, tmp$MSE_Pred))
                }
              }
            } else {
              for(i in seq(along=obj@spp)){
                if(isFALSE(obj@monotone[i]) & isTRUE(mono)){
                  suppD <- attr(obj@monotone, "D005")[i] * SKPar[[ SKspp[i] ]]$q001
                  suppH <- 0.01*obj@Ht[i]
                } else {
                  suppD <- suppH <- numeric(0)
                }
                tmp <- TapeR::E_DHx_HmDm_HT.f(Hx=Hx[ use[i] ],
                                              Dm = c(suppD, obj@Dm[[i]]),
                                              Hm = c(suppH, obj@Hm[[i]]),
                                              mHt = obj@Ht[i],
                                              sHt = obj@sHt[i],
                                              par.lme = SKPar[[ SKspp[i] ]],
                                              Rfn = Rfn)
                if(interval=="none"){
                  brk <- bark(Ba = rep(obj@spp[i], length(Hx[ use[i] ])),
                              Dm = tmp$DHx,
                              relH = Hx[ use[i] ]/obj@Ht[i])
                  res[i, ] <- tmp$DHx - brk
                } else if(interval=="confidence"){
                  brk <- bark(Ba = rep(obj@spp[i], length(Hx[ use[i] ]) * 3),
                              Dm = as.vector(tmp$CI_Mean),
                              relH = rep(Hx[ use[i] ]/obj@Ht[i], 3))
                  res <- rbind(res, cbind(tree=i,
                                          Hx[ use[i] ],
                                          tmp$CI_Mean - matrix(brk, ncol = 3, byrow = F)))
                } else if(interval=="prediction"){
                  brk <- bark(Ba = rep(obj@spp[i], length(Hx[ use[i] ]) * 3),
                              Dm = as.vector(tmp$CI_Pred),
                              relH = rep(Hx[ use[i] ]/obj@Ht[i], 3))
                  res <- rbind(res, cbind(tree=i,
                                          Hx[ use[i] ],
                                          tmp$CI_Pred - matrix(brk, ncol = 3, byrow = F)))
                } else if(interval=="MSE"){
                  brk <- bark(Ba = rep(obj@spp[i], length(Hx[ use[i] ])),
                              Dm = tmp$DHx,
                              relH = Hx[ use[i] ]/obj@Ht[i])
                  res <- rbind(res, cbind(tree=i, Hx[ use[i] ], tmp$DHx - brk,
                                          tmp$MSE_Mean, tmp$MSE_Pred))
                }
              }
            }
            colnames(res) <- colnms
            return(res[,,drop=TRUE])
          })
