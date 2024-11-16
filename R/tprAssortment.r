#' @title Functions to calculate assortments for given tree
#' @description Function calculates assortments for given tree according to
#' assortment specification
#' @param obj an object of class 'tprTrees'
#' @param pars parameters to specify assortments, see \code{\link{parSort}}
#' @param mono logical, defaults to true. If calibrated taper curve is
#' non-monotonic at stem base, a support diameter is added.
#' @param Rfn Rfn setting for residuals error matrix, defaults to
#' \code{list(fn="sig2")}, see \code{\link[TapeR]{resVar}}.
#' @return a data.frame with columns \code{tree}: tree identifier,
#' \code{sort}: assortment name, \code{height}: beginning of assortment along
#' trunk, \code{length}: length of assortment, \code{mdm}: mid-diameter of
#' assortment, \code{zdm}: top-diameter of assortment and \code{vol}: volume.
#' @import TapeR
#' @export
#' @examples
#' ## conifer wood
#' obj <- tprTrees(spp=c(1, 8),
#'                 Dm=list(30, 40),
#'                 Hm=list(1.3, 1.3),
#'                 Ht=c(30, 40))
#' tprAssortment(obj)
#' pars <- parSort(stH=0.2, Lxh=c(1, 1.5), fixN=2, fixL=4)
#' (ass <- tprAssortment(obj, pars))
#' plot(obj, assort = ass)
#'
#' ## deciduous wood
#' obj <- tprTrees(spp=c(15),
#'                 Dm=list(40),
#'                 Hm=list(1.3),
#'                 Ht=c(40))
#' tprAssortment(obj)
#' pars <- parSort(n=length(obj), Lxh=c(1), Hsh=10, Az=10)
#' ass <- tprAssortment(obj, pars)
#' plot(obj, assort=ass)
#'
setGeneric("tprAssortment",
           function(obj, pars=NULL, mono=TRUE, Rfn=NULL) standardGeneric("tprAssortment"))

#' @describeIn tprAssortment method for class 'tprTrees'
setMethod("tprAssortment", signature = "tprTrees",
          function(obj, pars=NULL, mono=TRUE, Rfn=NULL){
            if(is.null(pars)) pars=parSort(n=length(obj)) # default values

            if(!is.null(Rfn)){
              oldRfn <- options()$TapeS_Rfn
              options("TapeS_Rfn" = Rfn)
              on.exit(options("TapeS_Rfn" = oldRfn))
            }
            if(!identical(mono, oldmono <- options()$TapeS_mono)){
              options("TapeS_mono" = mono)
              on.exit(options("TapeS_mono" <- oldmono))
            }

            ## calculate total (physical volume)
            vol <- Vfm(obj) # 2m-segment-volume
            dbh <- Dbh(obj) # getting estimated Dbh

            ## check and update parameters and dependent variables
            pars@stH <- ifelse(pars@stH==0, 0.01*obj@Ht, pars@stH)
            obj@Ht <- ifelse(pars@Hkz==1, obj@Ht+2,
                             ifelse(pars@Hkz==2,
                                    ifelse(dbh<30, dbh, 30 + (dbh-30) * 0.3),
                                    obj@Ht))
            pars@Hsh <- ifelse(pars@Hsh==1 & obj@spp >= 15, 0.7*obj@Ht, pars@Hsh)
            pars@Hsh <- ifelse(pars@Hsh==2, 5, pars@Hsh) # TODO: Ndh maxH Sth und Ab = 5m
            pars@Hsh <- ifelse(pars@Hsh==3, pars@stH, pars@Hsh) # TODO: Ndh maxH Sth und Ab = 0.1m

            tmpAz <- Az(obj@spp, dbh)
            pars@Az <- ifelse(pars@Az==0, tmpAz, pars@Az)
            HAz <- tprHeight(obj, Dx = pars@Az, bark=TRUE, cp = FALSE)
            HAz <- ifelse(pars@Hsh==4 & HAz > 0.7*obj@Ht, 0.7*obj@Ht, HAz) # dead or broken stem

            pars@Zsh <- ifelse(pars@Zsh==0, pars@Az, pars@Zsh) # TODO: AZ mit bark vs. Zsh ohne bark?
            Hsh <- tprHeight(obj, Dx = pars@Zsh, bark=FALSE, cp = FALSE)
            Hsh <- ifelse(pars@Hsh != 0 & pars@Hsh < Hsh, pars@Hsh, Hsh)
            pars@Zab <- ifelse(pars@Zab==0, 14, pars@Zab)
            Hab <- tprHeight(obj, Dx = pars@Zab, bark=FALSE, cp = FALSE)
            Hab <- ifelse(pars@Hsh != 0 & pars@Hsh < Hab, pars@Hsh, Hab) # maximum height as given via pars@Hsh

            pars@trL <- ifelse(pars@trL == 0, 19, pars@trL)
            pars@fixL <- ifelse(pars@fixL > pars@trL, pars@trL, pars@fixL)# fixL no longer than transport length
            pars@fixZ <- ifelse(pars@fixZ==0, pars@Zsh, pars@fixZ)
            Hfix <- tprHeight(obj, Dx = pars@fixZ, bark = FALSE, cp = FALSE)
            Hfix <- ifelse(Hfix > HAz, HAz, Hfix)
            pars@fixA <- ifelse(pars@fixA>0, pars@fixA/100, pars@fixA) # transform from cm to m
            pars@fixA <- ifelse(pars@fixR==0 & pars@fixA==0, 0.01*pars@fixL, pars@fixA) # default 1% of fixL, given in m
            pars@fixA <- ifelse(pars@fixR*0.01*pars@fixL > pars@fixA,
                                pars@fixR*0.01*pars@fixL, pars@fixA)

            ## calculate assortments
            # xh = unusable wood at stem foot / X-Holz
            # fx = fix length at stem foot / Fixlängen
            # sh = stem wood / Stammholz
            # ab = top segment / Abschnitt
            # ih = industrial wood / Industrieholz
            # nvdh = unusable wood at top / nicht-verwertbares Derbholz
            assort <- list()
            for(i in seq(along=obj@spp)){
              # i <- 1
              acth <- pars@stH[i] # actual height

              ## calculate xh
              lxh <- pars@Lxh[i]
              if(lxh>0){
                hdx <- acth + 0.5*lxh
                dxhmr <- tprDiameter(obj[i], Hx = hdx, bark = TRUE)
                dxhmr <- dxhmr - ifelse(dxhmr >= 20, 0.75, 0.5)
                mdmxh <- dxhmr - bark(obj@spp[i], Dm = dxhmr, relH = hdx/obj@Ht[i])
                volxh <- pi * (mdmxh/200)^2 * lxh
                acthxh <- acth
                acth <- acth + lxh
              } else {
                acthxh <- 0
                volxh <- 0
                mdmxh <- 0
              }

              ## calculate fix length assortment
              if(pars@fixN[i] >= 1){
                actl <- Hfix[i] - acth # actual available length
                nfix <- as.integer(min(c(actl %/% pars@fixL[i], pars@fixN[i])))
                lfix <- pars@fixL[i] + pars@fixA[i] # gross length
                hfx <- acth + lfix/2 + (0:(nfix-1))*lfix # height of measurement
                dfxmr <- tprDiameterCpp(obj[i], Hx=hfx, bark = TRUE)
                dfxmr <- dfxmr - ifelse(dfxmr >= 20, 0.75, 0.5)
                mdmfx <- dfxmr - bark(rep(obj@spp[i], length(dfxmr)), Dm = dfxmr, relH = hfx/obj@Ht[i])
                lfx <- pars@fixL[i]
                volfx <- pi * (mdmfx/200)^2 * lfx
                acthfx <- acth + (0:(nfix-1) * lfix)
                acth <- acth + nfix * lfix
              } else {
                volfx <- 0
                mdmfx <- 0
                lfx <- 0
                acthfx <- 0
              }

              ## calculate sh
              if(acth + 3 < Hsh[i]){ # Sth should be at least 3m
                lsh <- floor((Hsh[i] - acth) * 10) / 10 # round to 10cm
                lsh <- ifelse(lsh > pars@trL[i], pars@trL[i], lsh)
                dshmr <- tprDiameterCpp(obj[i], Hx = acth + lsh/2, bark=TRUE)
                dshmr <- dshmr - ifelse(dshmr >= 20, 0.75, 0.5)
                mdmsh <- dshmr - bark(obj@spp[i], Dm = dshmr, relH = (acth + lsh/2)/obj@Ht[i])
                volsh <- pi * (mdmsh/200)^2 * lsh
                acthsh <- acth
                acth <- acth + lsh * 1.01 # TODO: variable add-on
              } else {
                volsh <- 0
                mdmsh <- 0
                lsh <- 0
                acthsh <- 0
              }

              ## calculate ab
              if(acth + 3 < Hab[i]){
                lab <- floor((Hab[i] - acth) * 10) / 10 # round to 10cm
                lab <- ifelse(lab > pars@trL[i], pars@trL[i], lab)
                dabmr <- tprDiameter(obj[i], Hx = acth + lab/2, bark=TRUE)
                dabmr <- dabmr - ifelse(dabmr >= 20, 0.75, 0.5)
                mdmab <- dabmr - bark(obj@spp[i], Dm = dabmr, relH = (acth + lab/2)/obj@Ht[i])
                volab <- pi * (mdmab/200)^2 * lab
                acthab <- acth
                acth <- acth + lab * 1.01 # TODO: variable add-on
              } else {
                volab <- 0
                mdmab <- 0
                lab <- 0
                acthab <- 0
              }

              ## calculate ih + nvdh
              if(obj@spp[i]<15){
                # conifers
                if((acth + 1) <= HAz[i]){ # at least 1m for industrial wood
                  lih <- floor(HAz[i] - acth) # round to full meter
                  lih <- ifelse(lih > pars@trL[i], pars@trL[i], lih)
                  if(pars@LIh[i] >=1){
                    ## fixed length without add on # TODO
                    nih <- lih %/% pars@LIh[i]
                    lih <- pars@LIh[i]
                    hih <- acth + lih/2 + (0:(nih-1)) * lih
                    dihmr <- tprDiameter(obj[i], Hx = hih, bark=TRUE)
                    dihmr <- dihmr - ifelse(dihmr >= 20, 0.75, 0.5)
                    mdmih <- dihmr - bark(rep(obj@spp[i], length(hih)), Dm = dihmr, relH = hih/obj@Ht[i])
                    volih <- pi * (mdmih/200)^2 * lih
                    acthih <- acth + (0:(nih-1)) * lih
                    acth <- acth + nih * lih * 1.0 # TODO: Ih variable add-on; in BDAT NO add-on!
                  } else {
                    # unrestricted length if pars@LIh == 0
                    dihmr <- tprDiameter(obj[i], Hx = acth + lih/2, bark=TRUE)
                    dihmr <- dihmr - ifelse(dihmr >= 20, 0.75, 0.5)
                    mdmih <- dihmr - bark(obj@spp[i], Dm = dihmr, relH = (acth + lab/2)/obj@Ht[i])
                    volih <- pi * (mdmih/200)^2 * lih
                    acthih <- acth
                    acth <- acth + lih * 1.0 # TODO: Ih variable add-on; in BDAT NO add-on!
                  }
                } else {
                  volih <- 0
                  mdmih <- 0
                  lih <- 0
                  acthih <- 0
                }
                ## calculate nvdh
                lnvd <- tprHeight(obj[i], Dx = 7, bark=TRUE, cp = FALSE) - acth
                dnvdmr <- tprDiameter(obj[i], Hx = acth + lnvd/2, bark = TRUE)
                dnvdmr <- dnvdmr - ifelse(dnvdmr >= 20, 0.75, 0.5)
                mdmnvd <- dnvdmr - bark(obj@spp[i], Dm = dnvdmr, relH = (acth + lnvd/2)/obj@Ht[i])
                volnvd <- pi * (mdmnvd/200)^2 * lnvd
                acthnvd <- acth

              } else {

                # deciduous wood
                HDhG <- tprHeight(obj[i], Dx = 7, bark=TRUE)
                if(TRUE){
                  acthih <- acth
                  volhshmr <- tprVolume(obj[i], AB = list(A=0, B=acth, sl=2),
                                        iAB = "H", bark = TRUE)
                  volkromr <- vol[i] - volhshmr # Kronenderbholzvolumen mit bark
                  ## schätzung nvdh
                  volu <- volur <- 0
                  if(pars@Az[i] > 8){ # like in BDAT
                    if(pars@Az[i] > 40) pars@Az[i] <- 40
                    if(dbh[i] > 60) d13 <- 60 else d13 <- dbh[i]
                    iAz <- floor(pars@Az[i] - 7) # Aufarbeitungszopf"klasse"
                    iD13 <- floor((d13-8)*0.5) + 1 # Durchmesserklasse
                    abwD <- d13-(iD13*2 + 6) # Abweichung zwischen Durchmesser und Durchmesser-Klasse [0; 2[
                    W1 <- fnUnvd(BaMap(obj@spp[i], 4), iD13, iAz)
                    if(iD13>26) W0 <- W1 else W0 <- fnUnvd(BaMap(obj@spp[i], 4), iD13+1, iAz)
                    unvdpr <- (W1-W0)*(2-abwD)/2 + W0
                    faktor <- vol[i] / volkromr * unvdpr / 100
                    if(faktor > 1) faktor <- 1
                    volur <- volkromr * faktor
                  }
                  volih <- 0
                  volihr <- volkromr - volur
                  ## schätzung mittl. Astdurchmesser
                  ka <- acth # Höhe Stammholzzopf = aktuelle Höhe (?!)
                  if(ka < 0.1 * obj@Ht[i]) ka <- 0.15 * obj@Ht[i] # unklar
                  if(BaMap(obj@spp[i], 5) == 1){ # Buchen-Modell
                    mADm <- 8.866 - .2721*ka + .2035*dbh[1] + .3546*obj@Ht[i]/ka + 1.4139*ka/dbh[i]
                  } else { # Eichen-Modell
                    mADm <- 6.473 - .119*ka + .1578*dbh[i] + 1.4236*obj@Ht[i]/ka
                  }
                  ## Ableitung Rindenanteil des Industrieholzes
                  if(abs(volihr)>0.0001){
                    Dacth <- tprDiameter(obj[i], Hx = acth)

                    if(HAz[i] < HDhG){
                      mADmih <- mADm + ((volkromr-volihr)/volkromr) * (Dacth-mADm)
                    } else {
                      mADmih <- mADm
                    }

                    HmADmih <- tprHeight(obj[i], Dx=mADmih, bark=TRUE)
                    Ri <- bark(obj@spp[i], mADmih, HmADmih / obj@Ht[i])

                    # Anteil Volumen ohne bark an Volumen mit bark
                    Riproz <- (mADmih)**2/(mADmih + Ri)**2

                    volih <- volihr * Riproz # Volumen ohne bark
                    mdmih <- 0 # mADmih # mittlere Astdurchmesser todo: alt: mdmih=0
                    lih <- 0 # 4 * volih / (pi * (mdmih/100)^2) # virtual length
                  } else {
                    volih <- 0
                    mdmih <- 0
                    lih <- 0
                  }

                } else {
                  volih <- 0
                  mdmih <- 0
                  lih <- 0
                }
                ## calculate nvdh
                acthnvd <- 0
                if(HAz[i] < HDhG){
                  mADmnvd <- mADm - ((volkromr-volur)/volkromr) * (mADm-7)
                  HmADmnvd <- tprHeight(obj[i], Dx=mADmnvd, bark=TRUE)
                  Ri <- bark(obj@spp[i], mADmnvd, HmADmnvd / obj@Ht[i])
                  Riproz <- mADmnvd**2 / (mADmnvd+Ri)**2
                  volnvd <- volur * Riproz
                  mdmnvd <- 0 # mADmnvd
                  lnvd <- 0 # 4 * volnvd / (pi * (mdmnvd/100)^2) # virtual length
                } else {
                  volnvd <- 0
                  mdmnvd <- 0
                  lnvd <- 0
                }

              }

              ## calculate EV
              volev <- vol[i] - (volxh + sum(volfx) + volsh + volab + sum(volih) + volnvd)

              ## build resulting data.frame holding all relevant assortment information
              sortfix <- paste0("fix", formatC(seq(along=volfx), width = 2, format = "d", flag = "0"))
              if(length(volih)> 1 | obj@spp[i]>14){
                sortih <- "ih"
              } else {
                sortih <- paste0("ih", formatC(seq(along=volih), width = 2, format = "d", flag = "0"))
              }

              df <- data.frame(tree=i,
                               sort=c("hx", sortfix, "sth", "ab", sortih, "nvd"),
                               height=c(acthxh, acthfx, acthsh, acthab, acthih, acthnvd),
                               length=c(lxh, rep(lfx, length(volfx)), lsh, lab, rep(lih, length(volih)), lnvd),
                               mdm=c(mdmxh, mdmfx, mdmsh, mdmab, mdmih, mdmnvd),
                               zdm=NA,
                               vol=c(volxh, volfx, volsh, volab, volih, volnvd))
              df$zdm <- tprDiameterCpp(obj[i], Hx = df$height + df$length, bark = TRUE)
              df$zdm <- df$zdm - ifelse(df$zdm >= 20, 0.75, 0.5)
              df$zdm <- df$zdm - bark(rep(obj@spp[i], nrow(df)), Dm = df$zdm,
                                       relH = (df$height + df$length)/obj@Ht[i])
              df$zdm <- ifelse(df$mdm==0, 0, df$zdm)
              df <- df[df$vol > 0.001, ]
              assort[[i]] <- df
            }
            res <- do.call(rbind, assort)

            return(res)
          })

