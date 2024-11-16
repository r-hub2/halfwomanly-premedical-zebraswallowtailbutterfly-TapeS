#' @title coerce different data sources into class 'tprTrees'
#' @description coercion functions to make NFI, segment and 'BDAT' data
#' available as 'tprTrees' objects
#' @details The coerced data is automatically checked for validity by the class
#' constructor. For available species codes see \code{\link{tprSpeciesCode}}.
#' @seealso \code{\link{tprTrees-class}}, \code{\link{tprTrees}},
#' \code{\link{tprSpeciesCode}}
#' @name data_coercion
#' @return an object of class 'tprTrees', see \code{\link{tprTrees-class}}
NULL

#' @param nfi data.frame with tree measurements as provided by german NFI
#' @param mapping mapping of column names
#' @details When coercing NFI (National forest inventory, BWI) data, one need
#' to provide the columns \code{BaTpr} (species code), \code{Bhd} (Dbh, [mm]),
#' \code{D03, [mm]} (diameter in 30\% of tree height) and \code{Hoehe} (tree
#' height, [dm]). Optionally, one can provide \code{H1} (measurement height of
#' \code{Bhd}, [dm]), \code{H2} (measurement height of \code{D03}, [dm]) as well
#' as \code{sHt} (measurement error of tree height, i.e. standard deviation [dm]);
#' otherwise these are assumed to be 1.3m, 30\% of tree height and 0 (zero),
#' respectively.
#'
#' Additionally, the NFI database stores diameter as [mm] and height as [dm]; it
#' is *not* necessary to transform to [cm] and [m], as the function does this.
#' Equally, \code{sHt} [dm] is transformed to \code{sHt} [m].
#'
#' Keep in mind that species codes of NFI are different from the taper models
#' for historical reasons (c.f. BDAT). Use the NFI table  ('x_Ba') to map
#' species codes beforehand (see examples).
#'
#' @describeIn data_coercion coercion of German NFI data
#' @export
#' @examples
#' # NFI data usually stored as integer and units: diameter=[mm] and height=[dm]
#' nfi <- data.frame(BaTpr=1L, Bhd=300L, D03=270L, Hoehe=250L)
#' tpr <- nfi_as_tprtrees(nfi)
#' tpr
#' tpr@sHt # defaults to 0
#'
#' # one can provide measurement heights explicitly
#' nfi <- data.frame(spp=1, Bhd=300, H1=12, D03=270, H=250)
#' nfi_as_tprtrees(nfi, mapping=c(spp="BaTpr", H="Hoehe"))
#'
#' # measurement error in height
#' nfi <- data.frame(BaTpr=1L, Bhd=300L, D03=270L, Hoehe=250L, sHt=15)
#' tpr <- nfi_as_tprtrees(nfi)
#' tpr@sHt

nfi_as_tprtrees <- function(nfi, mapping=NULL){

  if (is.data.frame(nfi)) {

    ## required columns
    reqcols <- c("BaTpr", "Bhd", "D03", "Hoehe", "H1", "H2")

    ## adjust data.frame to contain all necessary variables
    if (!("D03" %in% colnames(nfi)) & !("D03" %in% mapping))
      nfi$D03 <- 0
    if (!("H1" %in% colnames(nfi)) & !("H1" %in% mapping))
      nfi$H1 <- 13 # [dm]
    if (!("H2" %in% colnames(nfi)) & !("H2" %in% mapping))
      nfi$H2 <- NA
    if (!("sHt" %in% colnames(nfi)) & !("sHt" %in% mapping))
      nfi$sHt <- 0

    if (is.null(mapping)) {
      if (!all(reqcols %in% colnames(nfi))) {
        stop("either provide 'mapping' or name parameter 'nfi' appropriately!")
      }
    } else {
      if (!(is.character(mapping) && length(mapping) > 0)) {
        stop("'mapping' is not of type character or is not of appropriate length!")
      }
      for (i in seq(along = mapping)) {
        n <- which(colnames(nfi) == names(mapping)[i])
        colnames(nfi)[n] <- mapping[i]
      }
      if (!all(reqcols %in% colnames(nfi))) {
        stop("provided name mapping not sufficient!")
      }
    }
    # case, where H2 is not given, assume H2 belongs to D03
    if(all(is.na(nfi$H2))){
      nfi$H2 <- 0.3 * nfi$Hoehe
    }
  } else {
    stop("'nfi' not of class 'data.frame'")
  }

  return(tprTrees(spp = nfi$BaTpr,
                  Dm = lapply(seq(along=nfi$BaTpr), function(a){
                    d1 <- nfi$Bhd[a] / 10 # mm to cm
                    d2 <- nfi$D03[a] / 10 # mm to cm
                    c(d1[d1>0], d2[d2>0])}),
                  Hm = lapply(seq(along=nfi$BaTpr), function(a){
                    h1 <- nfi$H1[a] / 10 # dm to m
                    h2 <- nfi$H2[a] / 10 # dm to m
                    c(h1[nfi$Bhd[a]>0], h2[nfi$D03[a]>0])}),
                  Ht = nfi$Hoehe / 10,
                  sHt = nfi$sHt / 10))
}


#' @param seg data.frame with measured tree segments, see details.
#' @param mapping mapping of column names
#' @details Sectional measurements provide more information about the trunk of
#' a tree and are usually stored in a different way. They exhibit an arbitrary
#' amount of diameter measurements which also might vary from tree to tree.
#' Hence, \code{seg_as_tprtrees} expects a data.frame with columns
#' \code{Id}, \code{BaTpr} (species code), \code{Dm} (diameter measured, [cm]),
#' \code{Hm} (height of \code{Dm}, [m]) and optionally \code{Ht} (height of
#' tree, [m]). Tree height \code{Ht} can be included to
#' \code{Dm}-\code{Hm}-pairs with \code{Dm} being zero
#' (e.g. \code{Dm}=0, \code{Hm}=25). If \code{Ht} is given, it gains priority.
#' @describeIn data_coercion coercion of segmented data to class 'tprTrees'
#' @export
#' @examples
#' ## coercing sectional measurements
#' data(DxHx.df, package = "TapeR")
#' DxHx.df$BaTpr <- 1 # Norway spruce
#' segtprtrees <- seg_as_tprtrees(DxHx.df, mapping=c(Dx="Dm", Hx="Hm"))
#'
#' ## extract tree height from Dm-Hm measurements if not given explicitly
#' DxHx.df$Ht <- NULL # remove height, as already included with Dm=0
#' segtprtrees <- seg_as_tprtrees(DxHx.df, mapping=c(Dx="Dm", Hx="Hm"))
#' segtprtrees

seg_as_tprtrees <- function(seg, mapping = NULL){
  if (is.data.frame(seg)) {

    ## required columns
    reqcols <- c("Id", "BaTpr", "Dm", "Hm")
    ntrees <- length(unique(seg$Id))

    if (is.null(mapping)) {
      if (!all(reqcols %in% colnames(seg))) {
        stop("either provide 'mapping' or name parameter 'seg' appropriately!")
      }
    } else {
      if (!(is.character(mapping) && length(mapping) > 0)) {
        stop("'mapping' is not of type character or is not of appropriate length!")
      }
      for (i in seq(along = mapping)) {
        n <- which(colnames(seg) == names(mapping)[i])
        colnames(seg)[n] <- mapping[i]
      }
      if (!all(reqcols %in% colnames(seg))) {
        stop("provided name mapping not sufficient!")
      }
    }
  } else {
    stop("'seg' not of class 'data.frame'")
  }

  ## extract tree height if not explicitly given
  if (!("Ht" %in% colnames(seg))){
    ht <- seg[seg$Dm==0, c("Id", "Hm")]
    suId <- sort(unique(seg$Id))
    if((ntrees != nrow(ht)) | (!identical(suId, sort(ht$Id)))){
      stop("no tree height deducable for each tree-Id: either too many or too little obs with Dm==0!")
    } else {
      ht <- ht[order(ht$Id), "Hm", drop=TRUE]
    }
  } else {
    ht <- unique(seg[order(seg$Id), c("Id", "Ht")])$Ht
  }
  seg <- seg[seg$Dm > 0, ]
  return(tprTrees(spp = unique(seg[order(seg$Id), c("Id", "BaTpr")])$BaTpr,
                  Dm = lapply(sort(unique(seg$Id)), function(a){
                    seg[seg$Id==a, "Dm", drop=TRUE]
                    }),
                  Hm = lapply(sort(unique(seg$Id)), function(a){
                    seg[seg$Id==a, "Hm", drop=TRUE]
                    }),
                  Ht = ht))
}

#' @param bdat data.frame holding data to process with rBDAT
#' @details coercing object of class 'datBDAT' from R-Package "rBDAT" into
#' class 'tprTrees'
#' @describeIn data_coercion coercion of bdat data
#' @export
#' @examples
#' if(require(rBDAT)){
#'   bdt <- buildTree(list(spp=1, D1=30, D2=28, H2=7, H=25))
#'   bdat_as_tprtrees(bdt)
#' }
#'
#'
bdat_as_tprtrees <- function(bdat){
  if("datBDAT" %in% class(bdat)){
    if(any(bdat$D2 <= 0)){
      if(requireNamespace("rBDAT", quietly = TRUE)){
        bdat$Hx <- 0.3*bdat$H
        bdat$D2 <- rBDAT::getDiameter(bdat)
        bdat$H2 <- 0.3*bdat$H
      } else {
        bdat$D2 <- 0
        message("D2/H2 only considered if D2 > 0 or with 'rBDAT' installed \n D2 set to 0")
      }
    }
    bdat$H1 <- ifelse(bdat$H1 == 0, 1.3, bdat$H1)
    return(tprTrees(spp = bdat$spp,
                    Dm = lapply(seq(along=bdat$spp), function(a){
                      d <- c(bdat$D1[a], bdat$D2[a])
                      d[d>0]
                    }),
                    Hm = lapply(seq(along=bdat$spp), function(a){
                      d <- c(bdat$D1[a], bdat$D2[a])
                      h <- c(bdat$H1[a], bdat$H2[a])
                      h[d>0]
                    }),
                    Ht = bdat$H))
  } else {
    stop("'bdat' not of appropriate class")
  }
}
