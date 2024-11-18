#' @title Functions to calculate stem volume from taper curve
#' @description Function calculates stem volume from taper curve for given
#' trees, depending definition of segment and on bark indicator. It is possible
#' to request confidence or prediction intervals.
#' @param obj object of class 'tprTrees'
#' @param AB list with heights or diameters A and B of section for which volume
#' over or under bark should be calculated. Additionally, add in \code{sl} for
#' the segment length over which the integral should be calculated. See details.
#' @param iAB character indicating how to interpret given A and B values. Either
#' "H" (the default), "Dob" (diameter over bark) or "Dub" (diameter under bark).
#' Could be of length one or two, depending on whether A and B are both height
#' or diameter variables or not. See examples.
#' @param bark should volume be returned including (\code{TRUE}) or excluding
#' bark (\code{FALSE})?
#' @param interval character to indicate whether and which type of interval is
#' required; one of \code{none}, \code{confidence} or \code{prediction}.
#' @param mono logical, defaults to true. If calibrated taper curve is
#' non-monotonic at stem base, a support diameter is added.
#' @param Rfn Rfn setting for residuals error matrix, defaults to
#' \code{list(fn="sig2")}, see \code{\link[TapeR]{resVar}}.
#' @details The function returns total solid wood w/ bark (i.e. from H=0 to
#' D=7cm) by default. Using \code{AB}, one can specify lower \code{A} and upper
#' \code{B} end of segments for which volume is required, w/ or w/o bark.
#'
#' \code{iAB} can be a vector of length two, indicating how to interpret A and
#' B. Hence, one can calculate volume between a given height and a given
#' diameter, either over or under bark. If of length one, it is assumed that the
#' indicator applies to both A and B.
#'
#' Defining interval \code{'confidence'} or \code{'prediction'} returns
#' lower (lwr) and upper (upr) interval bounds on confidence level
#' \eqn{\alpha} = \code{qt(0.025, ...)}. NB: The volume confidence bounds only
#' incorporate the uncertainty of diameter estimation at a pre-fixed position
#' (e.g. H=1.3m). If the position is given as diameter (e.g. \code{iAB="Dob"}),
#' the absolute height position is calculated using the *estimated* diameter,
#' hence, the uncertainty of the estimated absolute height is not (yet)
#' included. Neither is the uncertainty of the models for bark reduction.
#'
#' In contrast to the underlying R-package TapeR, which uses
#' \code{\link[TapeR]{E_VOL_AB_HmDm_HT.f}} for volume calculation, this function
#' calculates volume based on stem-section (default: 2m, see parameter
#' \code{AB}). Additionally, with that approach, bark reduction is easily
#' possible.
#' @return if \code{interval='none'} a vector else a matrix.
#' @seealso \code{\link[TapeR]{E_DHx_HmDm_HT.f}} for the underlying diameter
#' calculation.
#' @import TapeR
#' @export
#' @examples
#' obj <- simTrees() # default is: simulate 10 Norway spruce with mean dbh of 40
#' A <- 1
#' B <- 10
#' tprVolume(obj) # default is: coarse wood volume w/ bark
#' tprVolume(obj, AB = list(A=A, B=B, sl=2), iAB = "H", bark=FALSE)
#' tprVolume(obj, AB = list(A=A, B=B, sl=0.01), iAB = "H", bark=FALSE)
#' tprVolume(obj, AB = list(A=A, B=B, sl=0.01), iAB = "H", bark=TRUE)
#'
#' ## compare against integrated taper curve volume via package TapeR
#' ## TapeR integrates over the taper curve, while TapeS uses segments of length 'sl'
#' SKP <- TapeS:::SKPar
#' TapeR::E_VOL_AB_HmDm_HT.f(Hm=obj@Hm[[1]], Dm = obj@Dm[[1]], iDH = "H",
#'                           mHt = obj@Ht[1], sHt = 0, A = A, B = B,
#'                           par.lme=SKP[[1]])$E_VOL
#'
#' ## returning intervals
#' tprVolume(obj, interval="none")
#' tprVolume(obj, interval="confidence")
#' tprVolume(obj, interval="prediction")
#' tprVolume(obj, interval="prediction", bark=FALSE)
#' tprVolume(obj, interval="prediction", AB=list(A=0.1, B=5.1, sl=0.1), iAB="H")

setGeneric("tprVolume",
           function(obj, AB=NULL, iAB=NULL, bark=NULL, interval="none", mono=TRUE, Rfn=NULL)
             standardGeneric("tprVolume"))

#' @describeIn tprVolume method for class 'tprTrees'
setMethod("tprVolume", signature = "tprTrees",
          function(obj, AB=list(A=0, B=7, sl=2), iAB=c("h", "dob"),
                   bark=TRUE, interval="none", mono=TRUE, Rfn=NULL){

            ## resulting object
            if(interval %in% c("confidence", "prediction")){
              ncol <- 4
              nrow <- length(obj)
              colnms <- c("lwr", "EVol", "upr", "MSE")
              if(identical(interval, "confidence")){
                useLE <- "DHx" # LE=list element
                useMSE <- "KOV_Mean" #
              } else {
                useLE <- "DHx" # prediction intervals
                useMSE <- "KOV_Pred" #
              }
            } else if(interval == "none"){
              ncol <- 1
              nrow <- length(obj)
              colnms <- NULL
              useLE <- "DHx"
              useMSE <- "KOV_Mean" # FIX: which interval required? which variance to be used?
            } else {
              stop("'interval' must be one of 'none', 'confidence', 'prediction'")
            }
            res <- matrix(NA, ncol = ncol, nrow = nrow, dimnames = list(NULL, colnms))

            ## get species code for taper curve / Schaftkurve
            SKspp <- BaMap(obj@spp, 1)
            if(is.null(Rfn)) Rfn <- getOption("TapeS_Rfn")

            ## check iAB: which element not valid? wenv
            wenv <- iAB[which(!(tolower(iAB) %in% c("h", "dub", "dob")))]
            if(length(wenv) > 0){
              stop(paste0("iAB should be in 'h', 'dub' or 'dob', not: ", wenv))
            }

            ## (re)define A, B and sl
            ## transform A into height if necessary
            iA <- tolower(iAB[1])
            iB <- tolower(ifelse(length(iAB)>1, iAB[2], iAB[1]))
            ## assert that A and B are appropriate (type, length, value)
            for(idx in c("A", "B")){
              stopifnot(length(AB[[idx]]) == 1 | length(AB[[idx]] == length(obj@Ht)),
                        all(AB[[idx]] >= 0), is.numeric(AB[[idx]]))
              if(length(AB[[idx]]) != length(obj@Ht)){
                # then length(AB[[idx]])==1 & length(obj@Ht)>1 and vector should be extended
                AB[[idx]] <- rep(AB[[idx]], length(obj@Ht))
              }
            }
            ## re-calculate
            if(iA == "dob"){
              A <- sapply(seq(length(obj@spp)), function(a){
                if(isFALSE(obj@monotone[a]) & isTRUE(mono)){
                  suppD <- attr(obj@monotone, "D005")[a] * SKPar[[ SKspp[a] ]]$q001
                  suppH <- 0.01*obj@Ht[a]
                } else {
                  suppD <- suppH <- numeric(0)
                }
                E_HDx_HmDm_HT.f(Dx = AB[["A"]][a],
                                Dm = c(suppD, obj@Dm[[a]]),
                                Hm = c(suppH, obj@Hm[[a]]),
                                mHt = obj@Ht[a],
                                sHt = 0,
                                par.lme = SKPar[[ SKspp[a] ]],
                                Rfn = Rfn)
              })
            } else if(iA == "dub"){
              A <- sapply(seq(length(obj@spp)), function(a){
                if(isFALSE(obj@monotone[a]) & isTRUE(mono)){
                  suppD <- attr(obj@monotone, "D005")[a] * SKPar[[ SKspp[a] ]]$q001
                  suppH <- 0.01*obj@Ht[a]
                } else {
                  suppD <- suppH <- numeric(0)
                }
                E_HDxoR_HmDm_Ht.f(DxoR = AB[["A"]][a],
                                  Dm = c(suppD, obj@Dm[[a]]),
                                  Hm = c(suppH, obj@Hm[[a]]),
                                  mHt = obj@Ht[a],
                                  sHt = 0,
                                  par.lme = SKPar[[ SKspp[a] ]],
                                  Rfn = Rfn)
              })
            } else {
              A <- AB[["A"]]
            }

            if(iB == "dob"){
              B <- sapply(seq(length(obj@spp)), function(a){
                if(isFALSE(obj@monotone[a]) & isTRUE(mono)){
                  suppD <- attr(obj@monotone, "D005")[a] * SKPar[[ SKspp[a] ]]$q001
                  suppH <- 0.01*obj@Ht[a]
                } else {
                  suppD <- suppH <- numeric(0)
                }
                TapeR::E_HDx_HmDm_HT.f(Dx = AB[["B"]][a],
                                       Dm = c(suppD, obj@Dm[[a]]),
                                       Hm = c(suppH, obj@Hm[[a]]),
                                       mHt = obj@Ht[a],
                                       sHt = 0,
                                       par.lme = SKPar[[ SKspp[a] ]],
                                       Rfn = Rfn)
              })
            } else if(iB == "dub"){
              B <- sapply(seq(length(obj@spp)), function(a){
                if(isFALSE(obj@monotone[a]) & isTRUE(mono)){
                  suppD <- attr(obj@monotone, "D005")[a] * SKPar[[ SKspp[a] ]]$q001
                  suppH <- 0.01*obj@Ht[a]
                } else {
                  suppD <- suppH <- numeric(0)
                }
                E_HDxoR_HmDm_Ht.f(DxoR = AB[["B"]][a],
                                  Dm = c(suppD, obj@Dm[[a]]),
                                  Hm = c(suppH, obj@Hm[[a]]),
                                  mHt = obj@Ht[a],
                                  sHt = 0,
                                  par.lme = SKPar[[ SKspp[a] ]],
                                  Rfn = Rfn)
              })
            } else {
              B <- AB[["B"]]
            }
            ## extract segment length (sl)
            sl <- ifelse(is.null(AB[["sl"]]), 2, AB[["sl"]]) # defaults=2m sections

            ## do for all trees
            if(isFALSE(bark)){
              ## excluding bark / under bark
              for(i in seq(along=obj@spp)){
                ab <- slfun(A[i], B[i], sl=sl)
                if(isFALSE(obj@monotone[i]) & isTRUE(mono)){
                  suppD <- attr(obj@monotone, "D005")[i] * SKPar[[ SKspp[i] ]]$q001
                  suppH <- 0.01*obj@Ht[i]
                } else {
                  suppD <- suppH <- numeric(0)
                }
                dlist <- TapeR::E_DHx_HmDm_HT.f(Hx = ab$Hx,
                                                Dm = c(suppD, obj@Dm[[i]]),
                                                Hm = c(suppH, obj@Hm[[i]]),
                                                mHt = obj@Ht[i],
                                                sHt = obj@sHt[i],
                                                par.lme = SKPar[[ SKspp[i] ]],
                                                Rfn = Rfn)
                D <- dlist[[ useLE ]]
                kovD <- dlist[[ useMSE ]]
                varD <- dlist[[ useMSE ]]
                D <- D - as.vector(bark(rep(obj@spp[i], length(dlist$DHx)),
                                        Dm = dlist$DHx,
                                        relH = ab$Hx/obj@Ht[i]))
                #https://de.wikipedia.org/wiki/Erwartungswert
                #Erwartungswert_des_Produkts_von_nicht_stochastisch_unabhaengigen_Zufallsvariablen
                #Falls die Zufallsvariablen X und Y nicht stochastisch unabhängig sind,
                #gilt für deren Produkt:
                # E ( X Y ) = E ( X ) E ( Y ) + Cov ( X , Y )
                # Dabei ist Cov ( X , Y ) die Kovarianz zwischen X und Y
                # hier, mit X=D und Y=D,  also VAR(D).
                vol <- colSums(pi/4 * 1e-4 * (D^2 + diag(kovD)) * ab$L)

                if(identical(interval, "none")){
                  ## only volume required
                  (res[i,] <- vol)
                } else {
                  ## volume and variance required
                  kovVol <- calcVCOVsekVol(D, kovD, ab$L) # vcov of segments
                  # Variance of total volume
                  # https://de.wikipedia.org/wiki/Varianz_(Stochastik)#Summen_und_Produkte
                  # Var(X + Y) = VAR(X) + VAR(Y) + 2*COV(X, Y)
                  varVol <- sum(diag(kovVol)) + 2*sum(kovVol[upper.tri(kovVol)])
                  c_alpha = qt(p=0.025, df=SKPar[[ SKspp[i] ]]$dfRes, ncp=0, lower.tail = FALSE)
                  lwr <- vol - c_alpha*sqrt(varVol)
                  upr <- vol + c_alpha*sqrt(varVol)
                  (res[i,] <- c(lwr, vol, upr, varVol))
                }
              }

            } else {

              ## including bark / over bark
              for(i in seq(along=obj@spp)){
                ab <- slfun(A[i], B[i], sl=sl)
                if(isFALSE(obj@monotone[i]) & isTRUE(mono)){
                  suppD <- attr(obj@monotone, "D005")[i] * SKPar[[ SKspp[i] ]]$q001
                  suppH <- 0.01*obj@Ht[i]
                } else {
                  suppD <- suppH <- numeric(0)
                }
                dlist <- TapeR::E_DHx_HmDm_HT.f(Hx = ab$Hx,
                                                Dm = c(suppD, obj@Dm[[i]]),
                                                Hm = c(suppH, obj@Hm[[i]]),
                                                mHt = obj@Ht[i],
                                                sHt = obj@sHt[i],
                                                par.lme = SKPar[[ SKspp[i] ]],
                                                Rfn = Rfn)
                D <- dlist[[ useLE ]]
                kovD <- dlist[[ useMSE ]] # not rounded (as compared to MSE_Mean/MSE_Pred)
                vol <- colSums(pi/4 * 1e-4 * (D^2 + diag(kovD)) * ab$L)

                if(identical(interval, "none")){
                  ## only volume required
                  (res[i,] <- vol)
                } else {
                  ## volume and variance required
                  kovVol <- calcVCOVsekVol(D, kovD, ab$L) # vcov of segments
                  # Variance of total volume
                  # https://de.wikipedia.org/wiki/Varianz_(Stochastik)#Summen_und_Produkte
                  # Var(X + Y) = VAR(X) + VAR(Y) + 2*COV(X, Y)
                  varVol <- sum(diag(kovVol)) + 2*sum(kovVol[upper.tri(kovVol)])
                  c_alpha = qt(p=0.025, df=SKPar[[ SKspp[i] ]]$dfRes, ncp=0, lower.tail = FALSE)
                  lwr <- vol - c_alpha*sqrt(varVol)
                  upr <- vol + c_alpha*sqrt(varVol)
                  (res[i,] <- c(lwr, vol, upr, varVol))
                }
              }
            }

            return(res[,,drop=TRUE])
          })


slfun <- function(A, B, sl){
  # a <- A - A # set to 0
  # b <- B - A
  # a <- 0
  # b <- 5.5
  # sl <- 1
  (nsek <- floor((B-A) / sl))
  if(nsek >= 1){ # reguläre Sektionen
    (regsek <- A + (0:(nsek-1)) * sl) # reguläre Sektionen
    (messregsek <- regsek + sl/2) # Messstellen reguläre Sektionen
  } else {
    messregsek <- numeric(0)
  }
  if((A + nsek*sl) != B){
    (maxregsek <- A + nsek * sl)
    (lastsek <- B - maxregsek)
    (messlastsek <- B - lastsek/2)
  } else {
    (lastsek <- numeric(0))
    (messlastsek <- numeric(0))
  }
  (messall <- c(messregsek, messlastsek))
  (seklen <- c(rep(sl, length(messregsek)), lastsek))
  return(list(Hx=messall, L=seklen))
}

