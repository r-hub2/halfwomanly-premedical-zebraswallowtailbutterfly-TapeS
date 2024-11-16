#' @title Functions to calculate double bark thickness for given diameter at
#' height \code{Hx}
#' @description Funktion evaluates the double bark thickness models developed
#' by Altherr et al (1974-79).
#' @param obj object of class 'tprTrees'
#' @param Hx height for which double bark thickness is required
#' @param cp cartesian product, i.e. apply all \code{Hx} to all trees, defaults
#' to TRUE
#' @param Rfn Rfn setting for residuals error matrix, defaults to
#' \code{list(fn="sig2")}, see \code{\link[TapeR]{resVar}}.
#' @param mono logical, defaults to true. If calibrated taper curve is
#' non-monotonic at stem base, an support diameter is added.
#' @return double bark thickness [cm]
#' @references Altherr, E., P. Unfried, J. Hradetzky and V. Hradetzky (1974).
#' Statistische Rindenbeziehungen als Hilfsmittel zur Ausformung und Aufmessung
#' unentrindeten Stammholzes. Kiefer, Buche, Hainbuche, Esche und Roterle.
#' Freiburg i. Br., Forstl. Versuchs- u. Forschungsanst. Baden-Württenberg.
#'
#' Altherr, E., P. Unfried, J. Hradetzky and V. Hradetzky (1975). Statistische
#' Rindenbeziehungen als Hilfsmittel zur Ausformung und Aufmessung unentrindeten
#'  Stammholzes. Europäische Lärche, Japanische Lärche, Schwarzkiefer,
#'  Stieleiche, Traubeneiche, Roteiche, Bergahorn und Linde. Freiburg i. Br.,
#'  Forstl. Versuchs- u. Forschungsanst. Baden-Württenberg.
#'
#' Altherr, E., P. Unfried, J. Hradetzky and V. Hradetzky (1976). Statistische
#' Rindenbeziehungen als Hilfsmittel zur Ausformung und Aufmessung unentrindeten
#'  Stammholzes. Weymouthskiefer, Robinie, Bergulme, Birke, Marilandica-Pappel
#'  und Robusta-Pappel. Freiburg i. Br., Forstl. Versuchs- u. Forschungsanst.
#'  Baden-Württenberg.
#'
#' Altherr, E., P. Unfried, J. Hradetzky and V. Hradetzky (1978). Statistische
#' Rindenbeziehungen als Hilfsmittel zur Ausformung und Aufmessung unentrindeten
#'  Stammholzes. Fichte, Tanne, Douglasie und Sitka-Fichte. Freiburg i. Br.,
#'  Forstl. Versuchs- u. Forschungsanst. Baden-Württenberg.
#'
#' Altherr, E., P. Unfried, J. Hradetzky and V. Hradetzky (1979). Statistische
#' Rindenbeziehungen als Hilfsmittel zur Ausformung und Aufmessung unentrindeten
#'  Stammholzes. Neupotz-Pappel, Regenerata-Pappel, Kirsche, Spitzahorn,
#'  Feldahorn, Aspe, Weide, Flatterulme, Tulpenbaum u. Elsbeere.
#'  Freiburg i. Br., Forstl. Versuchs- u. Forschungsanst. Baden-Württenberg.
#' @import TapeR
#' @export
#' @examples
#' ## calculating bark thickness depends on diameter estimation and hence on the
#' ## assumed residual variance at calibration.
#' ## can be Rfn=list(fn="sig2") (default), i.e. EBLUP estimation from taper curve
#' ## or e.g. Rfn=list(fn="zero"), i.e. force taper curve through the given measurements
#' options("TapeS_Rfn") # "sig2", default in TapeS
#' tmp <- tprTrees()
#' Dm(tmp); Hm(tmp) # Dbh = D(Hx=1.3) = 30cm (measured)
#' Dbh(tmp) # estimated via EBLUP from taper curve
#' tprBark(tmp, Hx = c(1.3, 5)) # bark thickness corresponds to Dbh(tmp)
#' (d <- tprDiameter(tmp, Hx = c(1.3, 5), bark=TRUE)) ## predicted
#' bark(1, d[1], 1.3/30) # the same!
#' bark(1, d[2], 5/30) # the same!
#'
#' ## if using option TapeS_Rfn = list(fn="zero"), force taper curve through measurements
#' setTapeSoptions(Rfn = list(fn="zero"))
#' options()$TapeS_Rfn
#' tprBark(tmp, Hx = c(1.3, 5))
#' bark(1, 30, 1.3/30) # the same but different to above
#' bark(1, d[1], 1.3/30) # cf. above
#' bark(1, 28, 5/30) # the same but different to above
#' bark(1, d[2], 1.3/30) # cf. above

setGeneric("tprBark",
           function(obj, Hx, cp=TRUE, mono=TRUE, Rfn=NULL)
             standardGeneric("tprBark"))

#' @describeIn tprBark method for class 'tprTrees'
setMethod("tprBark", signature = "tprTrees",
          function(obj, Hx, cp=TRUE, mono=TRUE, Rfn=NULL){
            SKspp <- BaMap(obj@spp, 1)
            if(is.null(Rfn)) Rfn <- getOption("TapeS_Rfn")
            if(isTRUE(cp)){
              ## apply all Hx to all trees subsequently via Hx[TRUE] => all Hx
              use <- rep(TRUE, length(obj))
            } else {
              ## for each tree, use exactly one given Hx
              use <- seq(along=obj)
              # recycle Hx to be at least of obj-length
              Hx <- Hx * rep(1, max(length(obj)*length(Hx), length(Hx)))
            }
            ncol <- ifelse(all(use==TRUE), length(Hx), 1)
            res <- matrix(NA, nrow = length(obj@spp), ncol = ncol)
            for(i in seq(along=obj@spp)){
              if(isFALSE(obj@monotone[i]) & isTRUE(mono)){
                suppD <- attr(obj@monotone, "D005")[i] * SKPar[[ SKspp[i] ]]$q001
                suppH <- 0.01*obj@Ht[i]
              } else {
                suppD <- suppH <- numeric(0)
              }
              DHx <- TapeR::E_DHx_HmDm_HT.f(Hx=Hx[ use[i] ],
                                            Dm = c(suppD, obj@Dm[[i]]),
                                            Hm = c(suppH, obj@Hm[[i]]),
                                            mHt = obj@Ht[i],
                                            par.lme = SKPar[[ SKspp[i] ]],
                                            Rfn = Rfn)$DHx
              res[i, ] <- bark(rep(obj@spp[i], length(Hx[ use[i] ])), DHx,
                               relH = Hx[ use[i] ]/obj@Ht[i])
            }

            return(res[,, drop=TRUE])
          })
