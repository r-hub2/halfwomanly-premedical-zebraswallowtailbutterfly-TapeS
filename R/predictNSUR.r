#' estimate variance components for component biomass functions
#'
#' @param data data / predictors given for prediction by \code{\link{nsur}}
#'  incl. species code for component biomass function, see \code{\link{BaMap}}.
#' @param estBM estimated biomass components for which variance information is
#'  required, given as data.frame, possibly use \code{df[,,drop=FALSE]}
#' @param comp which components are required, see \code{\link{tprBiomass}}
#' @param interval either \code{none}, \code{confidence} or \code{prediction}
#' @param level Tolerance / confidence level, defaults 0.95
#' @param adjVarPar should the variance information be taken from stable models?
#'  defaults to TRUE
#' @param as.list Should the return value be a list or \code{rbind} to a
#'  data.frame? Defaults to TRUE.
#'
#' @return a data.frame with information on lower and upper bound of required
#'  interval as well as the (given) estimate and the respective mean squared
#'  error
#' @details
#' Estimates confidence and prediction intervals according to the methods
#' presented in Parresol (2001).
#'
#' In case, \code{adjVarPar = TRUE}, the models with instable variance estimates
#' like Silver fir, Scots pine, Maple and Ash are, firstly, fitted by Norway
#' spruce and European beech, respectively, and, secondly, adjusted to the
#' expected value of the species specific model by substracting the difference
#' to the first model. With that, more stable and imho more realistic confidence
#' and prediction intervals are given. True, this assumes comparability of the
#' variances between species.
#'
#'
#'
#' @examples
#' d1 <- seq(42, 56, 2)
#' h <- estHeight(d1, 1)
#' data <- data.frame(spp = 1:8, # from BaMap(1, 7)
#'                    dbh = d1,
#'                    ht = h,
#'                    sth = 0.01*h,
#'                    D03 = 0.8 * d1,
#'                    kl = 0.7 * h)
#' estBM <- as.data.frame( nsur(spp = data$spp,
#'                              dbh = data$dbh,
#'                              ht = data$ht,
#'                              sth = data$sth,
#'                              d03 = data$D03,
#'                              kl = data$kl) )
#' estBM$agb <- rowSums(estBM[, -which(colnames(estBM)=="id")])
#' comp = c("sw", "agb")
#' interval = "confidence"
#' level = 0.95
#' adjVarPar = TRUE
#' e1 <- TapeS:::NSURvar(data, estBM, comp, interval="confidence", level=0.95, adjVarPar = TRUE)
#' e2 <- TapeS:::NSURvar(data, estBM, comp, interval="confidence", level=0.95, adjVarPar = FALSE)
#'
#' \dontrun{
#' par(mfrow=c(1, 2))
#' plot(x = data$dbh, y = e1$agb_ECBM, main="adjusted Var-Parameter", pch=data$spp,
#'      ylim=c(0.5*min(e1$agb_ECBM), 1.2*max(e1$agb_ECBM)), las=1,
#'      ylab="estimated AGB", xlab = "DBH [cm]")
#' invisible(sapply(1:nrow(e1), function(a){
#'   # a <- 1
#'   # lines(x = rep(data$dbh[a], 2), y = c(e2$agb_lwr[a], e2$agb_upr[a]),
#'   #       col="blue", lwd=2)
#'   rect(xleft = data$dbh[a] - 0.1, xright = data$dbh[a] + 0.1,
#'        ybottom = e2$agb_lwr[a], ytop = e2$agb_upr[a], border = "blue")
#'   lines(x = rep(data$dbh[a], 2), y = c(e1$agb_lwr[a], e1$agb_upr[a]),
#'         col="red", lwd=2)
#' }))
#' legend("bottomright", legend=c("Fi", "Ta", "Kie", "Dgl", "Bu", "Ei", "BAh", "Es"), pch=1:8)
#'
#' ## prediction intervals
#' e1 <- TapeS:::NSURvar(data, estBM, comp, interval="prediction", level=0.95, adjVarPar = TRUE)
#' e2 <- TapeS:::NSURvar(data, estBM, comp, interval="prediction", level=0.95, adjVarPar = FALSE)
#'
#' plot(x = data$dbh, y = e1$agb_ECBM, main="adjusted Var-Parameter", pch=data$spp,
#'      ylim=c(0, 2*max(e1$agb_ECBM)), las=1,
#'      ylab="estimated AGB", xlab = "DBH [cm]")
#' invisible(sapply(1:nrow(e1), function(a){
#'   # a <- 1
#'   # lines(x = rep(data$dbh[a], 2), y = c(e2$agb_lwr[a], e2$agb_upr[a]),
#'   #       col="blue", lwd=2)
#'   rect(xleft = data$dbh[a] - 0.1, xright = data$dbh[a] + 0.1,
#'        ybottom = e2$agb_lwr[a], ytop = e2$agb_upr[a], border = "blue")
#'   lines(x = rep(data$dbh[a], 2), y = c(e1$agb_lwr[a], e1$agb_upr[a]),
#'         col="red", lwd=2)
#' }))
#' legend("topleft", legend=c("Fi", "Ta", "Kie", "Dgl", "Bu", "Ei", "BAh", "Es"), pch=1:8)
#'
#' ## one species, large diameter range
#' spp <- 1 # spruce
#' spp <- 5 # beech
#' spp <- 2 # silver fir
#' spp <- 8 # ash
#' d1 <- seq(7, 80, 2)
#' h <- estHeight(d1, spp)
#' data <- data.frame(spp = spp,
#'                    dbh = d1,
#'                    ht = h,
#'                    sth = 0.01*h,
#'                    D03 = 0.8 * d1,
#'                    kl = 0.7 * h)
#' estBM <- as.data.frame( nsur(spp = data$spp,
#'                              dbh = data$dbh,
#'                              ht = data$ht,
#'                              sth = data$sth,
#'                              d03 = data$D03,
#'                              kl = data$kl) )
#' estBM$agb <- rowSums(estBM[, -which(colnames(estBM)=="id")])
#' comp = c("sw", "agb")
#' interval = "confidence"
#' level = 0.95
#' adjVarPar = TRUE
#' e1 <- TapeS:::NSURvar(data, estBM, comp, interval="confidence", level=0.95, adjVarPar = TRUE)
#' e2 <- TapeS:::NSURvar(data, estBM, comp, interval="confidence", level=0.95, adjVarPar = FALSE)
#'
#' par(mfrow=c(1, 2))
#' plot(x = data$dbh, y = e1$agb_ECBM, main="adjusted Var-Parameter", pch=data$spp,
#'      ylim=c(0.5*min(e1$agb_ECBM), 1.2*max(e1$agb_ECBM)), las=1,
#'      ylab="estimated AGB", xlab = "DBH [cm]")
#' invisible(sapply(1:nrow(e1), function(a){
#'   # a <- 1
#'   # lines(x = rep(data$dbh[a], 2), y = c(e2$agb_lwr[a], e2$agb_upr[a]),
#'   #       col="blue", lwd=2)
#'   rect(xleft = data$dbh[a] - 0.1, xright = data$dbh[a] + 0.1,
#'        ybottom = e2$agb_lwr[a], ytop = e2$agb_upr[a], border = "blue")
#'   lines(x = rep(data$dbh[a], 2), y = c(e1$agb_lwr[a], e1$agb_upr[a]),
#'         col="red", lwd=2)
#' }))
#'
#'
#' ## prediction intervals
#' e1 <- TapeS:::NSURvar(data, estBM, comp, interval="prediction", level=0.95, adjVarPar = TRUE)
#' e2 <- TapeS:::NSURvar(data, estBM, comp, interval="prediction", level=0.95, adjVarPar = FALSE)
#'
#' plot(x = data$dbh, y = e1$agb_ECBM, main="adjusted Var-Parameter", pch=data$spp,
#'      ylim=c(0, 2*max(e1$agb_ECBM)), las=1,
#'      ylab="estimated biomass", xlab = "DBH [cm]")
#' invisible(sapply(1:nrow(e1), function(a){
#'   # a <- 1
#'   # lines(x = rep(data$dbh[a], 2), y = c(e2$agb_lwr[a], e2$agb_upr[a]),
#'   #       col="blue", lwd=2)
#'   rect(xleft = data$dbh[a] - 0.1, xright = data$dbh[a] + 0.1,
#'        ybottom = e2$agb_lwr[a], ytop = e2$agb_upr[a], border = "blue")
#'   lines(x = rep(data$dbh[a], 2), y = c(e1$agb_lwr[a], e1$agb_upr[a]),
#'         col="red", lwd=2)
#' }))
#' }
#'
#'

NSURvar <- function(data, estBM = NULL, comp = NULL, interval = "confidence",
                    level = 0.95, adjVarPar = TRUE, as.list=TRUE){

  if(is.null(estBM) & is.null(comp)){
    stop("either 'estBM' or 'comp' must be given!")
  } else if(!is.null(estBM) & is.null(comp)){
    comp <- colnames(estBM)
    comp <- check_Comp(comp)
  } else if(is.null(estBM) & !is.null(comp)){
    comp <- check_Comp(comp)
    estBM <- as.data.frame( nsur(spp=data$spp, dbh = data$dbh, ht = data$ht,
                                 sth = data$sth, d03 = data$D03, kl = data$kl))
    estBM$agb <- rowSums(estBM[, -which(colnames(estBM)=="id")])
  } else {
    # if both are given, check if all required comps are in estBM
    stopifnot("fn NSURvar: given 'comp' not in 'estBM'"=all(comp %in% colnames(estBM)))
  }

  if(!(interval[1] %in% c("confidence", "prediction"))){
    stop("given interval name unkown.")
  } else {
    interval <- interval[1]
  }

  # update data,
  # because in biomass function we need adjusted tree height (ht=ht-sth) and DH=dbh*ht
  # the order of the columns is given by the calling function tprBiomass()
  colnames(data) <- c("spp", "BHD", "Hoehe", "Stockhoehe", "D03", "KL")
  data$Hoehe <- data$Hoehe - data$Stockhoehe
  data$DH <- data$BHD * data$Hoehe
  data$no <- 1:nrow(data) # add row identifier to combine final data
  # head(data)

  uspp <- sort(unique(data$spp))
  bmsp <- c("fi", "ta", "kie", "dgl", "bu", "ei", "ba", "es")
  res <- list()
  for(u in uspp){
    # u <- 1
    # remove all variables with names of form 'axx' (=the parameters)
    rm(list=ls()[grepl("^a[1-9]{2}", x = ls())], inherits = FALSE)
    df <- data[data$spp == u,]
    est <- estBM[data$spp == u,, drop=FALSE]

    if(isTRUE(adjVarPar)){
      uadj <- ifelse(u %in% c(2, 4), 1, #replace for silver fir and scots pine by spruce
                     ifelse(u >= 7, 5, u)) # replace maple and ash by beech
      # BMM <- TapeS:::BMM
      obj <- BMM[[ bmsp[ uadj ] ]]
      # obj$qbar <- BMM[[ bmsp[ u ] ]]$qbar # still use the parameter estimates
      # obj$t <- BMM[[ bmsp[ u ] ]]$t # and uncertainty of current species
      # obj$formula <- BMM[[ bmsp[ u ] ]]$formula # and formula
      obj$n <- BMM[[ bmsp[ u ] ]]$n # n
      obj$k <- BMM[[ bmsp[ u ] ]]$k # k
    } else {
      obj <- BMM[[ bmsp[ u ] ]]
    }

    (EqNames <- obj$eqnlabel)
    if(length(EqNames)==8){
      names(EqNames) <- c("stw", "stb", "stwb", "sw", "sb", "swb", "fwb", "agb")
    } else if(length(EqNames)==9){
      names(EqNames) <- c("stw", "stb", "stwb", "sw", "sb", "swb", "fwb", "ndl", "agb")
    }
    selectedcomp <- EqNames[ match(comp, names(EqNames)) ]
    eq  <- which(EqNames %in% selectedcomp)
    eqnnms <- EqNames[eq]

    ## make parameters available
    ## take them from pooled results!!!
    for (i in 1:length(obj$qbar)) {
      name <- names(obj$qbar)[i]
      val <- obj$qbar[i]
      storage.mode(val) <- "double"
      assign(name, val)
    }

    if("prediction" %in% interval){
      ## calculate weights to be used in estimation of prediction interval
      Phi <- as.data.frame(matrix(0, nrow=nrow(df), ncol=length(eq)))
      colnames(Phi) <- eqnnms
      for(enm in eqnnms){
        # enm <- eqnnms[1]
        covar <- strsplit(as.character(obj$varModel[[enm]][[2]]), split = "|", fixed = T)
        covar <- ifelse(length(covar) > 1, covar[[2]], covar[[1]])
        Phi[, enm] <- df[, covar]^(obj$varPar[names(obj$varPar)==enm])
      }
    }

    ## prediction for each eqns
    resi <- list() # result for equations
    for(i in seq(along=eq)){
      # i <- 1
      ## given mean prediction
      (ECBM <- est[[names(EqNames[ eq[i] ])]])

      ## create some necessary variables
      (n <- obj$n[ eq[i] ]) # number of observation in each component
      (k <- obj$k[ eq[i] ]) # number of parameters
      (t <- qt(1-(1-level)/2, df = n-k))

      ## estimate confidence interval
      (fb <- as.matrix(attr(with(df, with(as.list(obj$qbar),
                                              eval(deriv(obj$formula[[ eq[i] ]],
                                                         names(obj$qbar))))), "gradient")))
      (estVar <- apply(fb, 1, function(a){t(a) %*% obj$t %*% a}))

      if("confidence" %in% interval){
        (lwr <- ECBM - (t * sqrt(estVar)))
        (MSE <- estVar) # estimated variance
        (upr <- ECBM + (t * sqrt(estVar)))
      }

      if("prediction" %in% interval){
        ## estimate prediction interval
        ## obs: obj$sigma_nsur holds the NSUR residual variance
        ## (in fitting it is called 'sysvar'; when pooling it becomes 'sigma_nsur')
        sii <- diag(obj$Sigma)[[ eqnnms[i] ]] ## covariance of the errors for equation i
        w <- Phi[, i] ## weight for equation i
        (MSE <- estVar + obj$sigma_nsur * sii * w)
        (lwr <- ECBM - (t * sqrt(MSE)))
        (upr <- ECBM + (t * sqrt(MSE)))
      }
      tmp <- data.frame(lwr, ECBM, upr, MSE)
      colnames(tmp) <- paste0(names(eqnnms)[i], "_", colnames(tmp))
      resi[[ i ]] <- tmp

    } #for loop equation

    if("ndl" %in% comp & any(1:4 %in% uspp) & u > 4){
      # there is no needle comp for deciduous tree species, but needed due to
      # at least one conifer in the data set
      resi <- cbind(resi, data.frame(ndl_lwr=NA,
                                     ndl_ECBM=NA,
                                     ndl_upr=NA,
                                     ndl_MSE=NA))
    }
    res[[u]] <- do.call(cbind, resi)
  } #for loop species

  res <- do.call(rbind, res)
  if(isTRUE(as.list)){
    ll <- list()
    for(cmp in comp){
      # cmp <- comp[1]
      ll[[cmp]] <- res[, grep(paste0("^", cmp), colnames(res))]
    }
    res <- ll
  }
  return(res)
}

