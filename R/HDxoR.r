#' @title Find Height of diameter under bark via uniroot
#' @description Functional equivalent to \code{\link[TapeR]{E_HDx_HmDm_HT.f}},
#' finding the height of a given diameter *without* bark, i.e. double bark
#' thickness needs to be added on top of given diameter to find appropriate
#' height.
#' @param DxoR Scalar. Diameter under bark for which to return height.
#' @param Hm Numeric vector of stem heights (m) along which diameter
#' measurements were taken for calibration. Can be of length 1. Must be of same
#' length as \code{Dm}.
#' @param Dm Numeric vector of diameter measurements (cm) taken for calibration.
#'  Can be of length 1. Must be of same length as \code{Hm}.
#' @param mHt Scalar. Tree height (m).
#' @param sHt Scalar. Standard deviation of stem height. Can be 0 if height was
#' measured without error.
#' @param par.lme List of taper model parameters obtained by
#' \code{\link[TapeR]{TapeR_FIT_LME.f}}, enhanced by the attribute 'spp', which
#' refers to the tree species used for double bark thickness
#' @param Rfn setting for residuals error matrix, defaults to \code{"sig2"}, see
#' details.
#' @param ... not currently used
#' @details finds height of given diameter via \code{uniroot}.
#' @return A scalar. Estimated height (m) given a diameter without bark.
#' @importFrom stats uniroot
#' @export
#' @examples
#' tmp <- tprTrees()
#' spp <- spp(tmp)
#' Hm <- Hm(tmp)
#' Dm <- Dm(tmp)
#' H <- Ht(tmp)
#' SKP <- TapeS:::SKPar
#' sppSK <- BaMap(spp, 1) # tree species for taper curve
#' ## diameter in 5m height
#' TapeR::E_DHx_HmDm_HT.f(c(5, 10), Hm, Dm, mHt=H, sHt = 0, par.lme = SKP[[sppSK]])$DHx
#' (D5m <- TapeR::E_DHx_HmDm_HT.f(c(5, 10), Hm, Dm, mHt=H, sHt = 0, par.lme = SKP[[sppSK]])$DHx)
#' ## bark thickness of diameter in 5m height
#' (RiD5m <- bark(c(1,1), Dm = D5m, relH = c(5, 10)/H))
#' ## find height of diameter without bark, which should be 5m
#' d5mub <- D5m - RiD5m
#' E_HDxoR_HmDm_Ht.f(DxoR = d5mub, Hm = Hm, Dm = Dm, mHt = H,
#'                   sHt = 0, par.lme = SKP[[sppSK]])

E_HDxoR_HmDm_Ht.f <- function (DxoR, Hm, Dm, mHt, sHt = 0, par.lme, Rfn=NULL, ...)
{
  if(is.null(Rfn)) Rfn <- getOption("TapeS_Rfn")
  xseq <- seq(0, mHt, length.out = 101)
  tpr <- tprTrees(spp = attr(par.lme, "spp"), Dm = Dm, Hm = Hm, Ht = mHt)

  Hx <- numeric(length(as.vector(DxoR)))
  for(k in seq(along=as.vector(DxoR))){
    d <- as.vector(DxoR)[k]
    Ddiff <- tprDiameter(tpr, Hx = xseq, bark=FALSE) - d
    xnull <- xseq[which(Ddiff == 0)]
    Dprod <- Ddiff[1:100] * Ddiff[2:101]
    for(i in which(Dprod < 0)){
      xnull <- c(xnull,
                 uniroot(HxoR_root.f, c(xseq[i], xseq[i + 1]), tol = 0.00001,
                         d, Hm, Dm, mHt, sHt, par.lme, Rfn)$root)
    }
    Hx[k] = max(xnull)
  }


  return(Hx)
}

#' @param Hx height at which taper curve is evaluated
#' @describeIn E_HDxoR_HmDm_Ht.f function to be searched
#'
HxoR_root.f <- function (Hx, DxoR, Hm, Dm, mHt, sHt, par.lme, Rfn, ...)
{
  DxmR <- TapeR::E_DHx_HmDm_HT.f(Hx, Hm, Dm, mHt, sHt = 0, par.lme, Rfn)$DHx
  RindeDxmR <- bark(Ba = attr(par.lme, "spp"), Dm = DxmR, relH = Hx/mHt)
  return(DxmR - (DxoR + RindeDxmR))
}
