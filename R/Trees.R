##########
# define class 'tprTrees'
##########

## check validity ####
check_trees <- function(object) {
  errors <- character()

  if (!is.integer(object@spp)) {
    msg <- "spp should of type 'integer'"
    errors <- c(errors, msg)
  }

  if (any(object@spp > 39) | any(object@spp < 1)) {
    msg <- "spp should be inside [1; 39]"
    errors <- c(errors, msg)
  }

  if(!is.numeric(object@Ht)) {
    msg <- "'Ht' should be of type 'numeric'"
    errors <- c(errors, msg)
  }

  if (any(object@Ht <= 0)) {
    msg <- "'Ht' should be positiv real"
    errors <- c(errors, msg)
  }

  if(!is.numeric(object@sHt)) {
    msg <- "'sHt' should be of type 'numeric'"
    errors <- c(errors, msg)
  }

  if (any(object@sHt < 0)) {
    msg <- "'sHt' should be positive real or zero"
    errors <- c(errors, msg)
  }

  if(!is.numeric(unlist(object@Dm))) {
    msg <- "'Dm' should be of type 'numeric'"
    errors <- c(errors, msg)
  }

  if(any(unlist(object@Dm) <= 0)) {
    msg <- "'Dm' should be positiv real"
    errors <- c(errors, msg)
  }

  if(!is.numeric(unlist(object@Hm))) {
    msg <- "'Hm' should be of type 'numeric'"
    errors <- c(errors, msg)
  }

  if(any(unlist(object@Hm) <= 0)) {
    msg <- "'Hm' should be positiv real"
    errors <- c(errors, msg)
  }

  if(any(sapply(1:length(object@spp), function(a){
    length(object@Hm[[a]]) != length(object@Dm[[a]])}) )){
    msg <- "'Dm' and 'Hm' should be of pairwise same length"
    errors <- c(errors, msg)
  }

  if(any(m <- sapply(1:length(object@spp), function(a){
    any(object@Ht[a] < object@Hm[[a]])}) )){
    msg <- paste0("'Hm' should be less than 'Ht' in each element")
    msg2 <- paste0("; check element i = ", paste0(which(m==TRUE), collapse = ", "))
    errors <- c(errors, paste0(msg, msg2))
  }

  if(!is.logical(object@monotone)){
    msg <- "'montone' should be logical"
    errors <- c(errors, msg)
  }

  if(length(errors) == 0){
    ## check monotonicity of taper curve
    object@monotone <- TRUE
  }

  if (length(errors) == 0){
    return(TRUE)
  } else {
    return(errors)
  }
}

## set class ####
#' @title An S4 class to represent one or multiple trees.
#' @description This class represents one or multiple trees by their biometric
#' characteristics.
#' @slot spp species code of trees
#' @slot Dm list of measured diameters
#' @slot Hm list of heights of measured diameters
#' @slot Ht total height of trees
#' @slot sHt standard deviation of total tree height, defaults to 0 for exact
#' height measurements without error
#' @slot monotone logical indicator about monotonicity of taper curve
#' @details blabla
#' @importFrom methods is new slot slotNames validObject
#' @examples
#' tprTrees() # initialise object by constructor
#' (tmp <- tprTrees(spp=c(1L,3L), Dm=list(c(30, 28), c(40, 38)),
#'                  Hm=list(c(1.3, 5), c(1.3, 5)), Ht=c(30, 40)))
#'
setClass("tprTrees",
         representation(spp = "integer",
                        Dm = "list",
                        Hm = "list",
                        Ht = "numeric",
                        sHt = "numeric",
                        monotone = "logical"),
         prototype(spp=1L, Dm=list(c(30, 28)), Hm=list(c(1.3, 5)), Ht=30, sHt=0,
                   monotone = FALSE),
         validity = check_trees)


## define subsetting ####
#' @title subsetting an object of class 'tprTree'
#' @description using indices i and j to subset
## #' @name [-method
#' @describeIn tprTrees subsetting for class 'tprTrees'
#' @aliases [,tprTrees-method
#' @aliases [,tprTrees,ANY,ANY,ANY-method
#' @keywords methods
#' @docType methods
#' @param x object from which to extract
#' @param i index i
#' @param j index j
#' @param ... not currently used
#' @param drop drop dimensions, defaults to FALSE
#' @return a part of the original object
#' @export
setMethod("[", signature(x = "tprTrees", i="ANY", j="ANY", drop="ANY"),
          function(x, i, j, ..., drop=FALSE){
            if(!is(x, "tprTrees"))
              stop("'x' needs to be of class 'tprTrees'")
            validObject(x)
            x@spp <- x@spp[i]
            x@Dm <- x@Dm[i]
            x@Hm <- x@Hm[i]
            x@Ht <- x@Ht[i]
            x@sHt <- x@sHt[i]
            attrbts1 <- attr(x@monotone, "Rfn")
            attrbts2 <- attr(x@monotone, "D005")
            x@monotone <- x@monotone[i]
            attr(x@monotone, "Rfn") <- attrbts1
            attr(x@monotone, "D005") <- attrbts2[i]
            if(missing(j)) return(x)
            if(is.numeric(j)){
              sN <- slotNames(x)
              x <- sapply(j, function(a) unlist(slot(x, sN[a])))
              # names(x) <- sN[j]
            } else {
              x <- t(sapply(j, function(a) unlist(slot(x, a))))
              # names(x) <- j
            }
            return(x)
          })

## define length function ####
#' @describeIn tprTrees length function for class 'tprTrees'
#' @param x object of class 'tprTrees'
setMethod("length", signature(x = "tprTrees"),
          function(x){
            if(!is(x, "tprTrees"))
              stop("'x' needs to be of class 'tprTree'")
            validObject(x)
            return(length(x@spp))
          })

## define 'show' method
#' @importFrom stats aggregate
#' @describeIn tprTrees length function for class 'tprTrees'
#' @param object object of class 'tprTrees'
setMethod("show", signature(object = "tprTrees"),
function(object){
    l <- aggregate(list(l=spp(object)), by=list(spp=spp(object)), FUN=length)$l
    cat("\nobject of class 'tprTrees'")
    cat("\n\ttotal number of trees:", sum(l))
    cat("\n\tspecies included: ", paste0(sort(unique(object@spp)),
                                             collapse = "\t"))
    cat("\n\t                  ",
        paste0(tprSpeciesCode(sort(unique(object@spp)), "kurz"),
               collapse = " \t"))
    cat("\n\tnumber of trees:  ", paste0(l, collapse = "\t"))
    dbh <- aggregate(list(d=Dbh(object)), by=list(spp=object@spp), FUN="mean")$d
    cat("\n\tmean dbh:         ",
        paste0(round(dbh, 1), collapse = "\t"))
    ht <- aggregate(list(h=Ht(object)), by=list(spp=object@spp), FUN="mean")$h
    cat("\n\tmean total height:",
        paste0(round(ht, 1), collapse = "\t"))
    nobs <- sapply(seq(along=object@spp), function(a){
      length(object@Dm[[a]])
      })
    nobs <- aggregate(list(n=nobs), by=list(spp=object@spp), FUN="mean")$n
    cat("\n\tmean obs/tree:    ",
        paste0(round(nobs, 1), collapse = "\t"))
    cat("\n\tmonotone taper curves: ", sum(object@monotone), "/", sum(l))
    cat("\n- ResVar-Matrix by fn '", attr(object@monotone, "Rfn")$fn, "'", sep="")
    if(sum(object@monotone) < sum(l)){
      cat("\n- consider using mono=TRUE (default) for taper curve evaluation ",
          "\n  if non-monotone taper curves are present ")
    }

})

## define 'summary' of object


## constructor ####
#' @title constructor for class tprTrees
#' @param spp species code, see \code{\link{tprSpeciesCode}}
#' @param Dm measurements of diameter along trunk
#' @param Hm height of measurements along trunk
#' @param Ht tree height
#' @param sHt standard deviation of stem height \code{Ht}. Can be 0 if height
#' was measured without error.
#' @param inv indicator (0-5) for inventory to assess taper form; numeric scalar
#' see \code{\link{FormTariff}}
#' @param Rfn function to populate residual variance matrix R
#' @param ... arguments to be passed to \code{initialize()}
#' @details constructor for a tprTrees object, includes a check on monotonicity
#' of the taper curve.
#' @return object of class \code{tprTrees}.
#' @examples
#' # just define a tree
#' tpr <- tprTrees(spp=1, Dm=30, Hm=1.3, Ht=27)
#' plot(tpr)
#' # define 2 trees with only dbh
#' tpr <- tprTrees(spp=c(1,3), Dm=c(30, 35), Hm=c(1.3, 1.3), Ht=c(27, 30))
#' plot(tpr)
#' # define 2 trees with several measurement
#' tpr <- tprTrees(spp=c(1,3), Dm=list(c(30, 28), c(35, 33, 31)),
#'                 Hm=list(c(1.3, 8), c(1.3, 5, 8)), Ht=c(27, 30))
#' plot(tpr)
#' # define 2 trees with only dbh and inventory indicator (form)
#' tpr <- tprTrees(spp=c(1,3), Dm=c(30, 35), Hm=c(1.3, 1.3), Ht=c(27, 30), inv=4)
#' plot(tpr)

#' @export
tprTrees <- function(spp=1L, Dm=list(c(30, 28)), Hm=list(c(1.3, 5)), Ht=30,
                     sHt=rep(0, length(Ht)), inv=NULL, Rfn=NULL, ...){

  if(!is.list(Dm) & !is.list(Hm) & length(Dm)==length(Hm)){
    if(length(Dm) == length(Ht)){
      ## same length, i.e. each element corresponds to one tree

      if(!is.null(inv) & all(round(Hm, 1) == 1.3)){
        ## additionally, form information available
        inv <- ifelse(inv < 0 | inv > 5, 0, inv)
        q03 <- FormTariff(spp, Dm, Ht, inv)
        # TODO: this is like Münchhausen
        # honestly, D005 is not exactly known, here we assume standard form
        # otherwise we would need to iterate to find the correct D005 w.r.t. Dbh and q03
        D005 <- tprDiameter(tprTrees(spp=spp, Dm=Dm, Hm=Hm, Ht=Ht),
                            Hx=0.05*Ht, cp=FALSE)
        Dm <- lapply(seq(length(Dm)), function(a){c(Dm[[a]], q03[a]*D005[[a]])})
        Hm <- lapply(seq(length(Hm)), function(a){c(Hm[[a]], 0.3*Ht[a])})

      } else {
        Dm <- as.list(Dm)
        Hm <- as.list(Hm)
      }

    } else {
      ## just coerce to list; testing follows in constructor
      Dm <- list(Dm)
      Hm <- list(Hm)
    }
  }

  x <- new("tprTrees",
           spp=as.integer(spp),
           Dm=Dm,
           Hm=Hm,
           Ht=Ht,
           sHt=sHt,
           monotone = rep(FALSE, length(spp)),
           ...)
  if(validObject(x)){
    x@monotone <- check_monotonicity(x, Rfn=Rfn)
  }
  return(x)
}

#' @title monotonicity check for taper curve
#' @param obj object of class 'tprTrees'
#' @param Rfn Rfn setting for residuals error matrix, defaults to
#' \code{list(fn="sig2")}, see \code{\link[TapeR]{resVar}}.
#' @details Taper curves are required to decrease monotonically. To avoid the
#' evaluation of non-monotone taper curves, a check is done through the
#' constructor function and an indicator (\code{monotone}) is set for each tree
#' stored inside the \code{tprTrees}-class. As the data has been check on validity
#' before this function is applied, we can use the tpr*-functions to evaluate
#' the taper curve and its monotonicity.
#' The check is done via comparison of the expected diameters along the trunk in
#' 1m-steps and its sorted (monotonically decreasing) version using
#' \code{\link{identical}}.
#' @return vector of logicals, same length as \code{spp}.
#' @export
check_monotonicity <- function(obj, Rfn=NULL){

  ## evaluate taper curve using cpp-function lmeSKEBLUP
  SKspp <- BaMap(obj@spp, 1)
  if(is.null(Rfn)) Rfn <- getOption("TapeS_Rfn")
  mono <- lapply(seq(along = obj@spp), function(i){
    Hx <- seq(0, 0.3, by = .01)
    rv <- TapeR::resVar(obj@Hm[[i]]/obj@Ht[i], fn=Rfn$fn,
                         sig2=SKPar[[ SKspp[i] ]]$sig2_eps, par=Rfn$par)
    d <- as.vector(lmeSKEBLUP(xm = obj@Hm[[i]]/obj@Ht[i],
                              ym = obj@Dm[[i]],
                              xp = Hx,
                              par = SKPar[[ SKspp[i] ]],
                              RV = rv)$yp)
    ## check monotonicity, but allow for diameter increase of max 1% of diameter
    TF <- identical(d, sort(d, decreasing = TRUE)) # | all(diff(d) / d[-1] <= .01)
    lst <- list(TF=TF
                # , maxD=max(d)
                # , HmaxD=Hx[which.max(d)]
                , D005=d[which(Hx==0.05)])
    return(lst)
  })
  res <- sapply(1:length(mono), function(a) mono[[a]]$TF)
  attr(res, "Rfn") <- options()$TapeS_Rfn
  attr(res, "D005") <- sapply(1:length(mono), function(a) mono[[a]]$D005)
  return(res)
}
