##########
# define class parSort
##########

## check validity ####
check_parSort <- function(object) {
  errors <- character()

  ## stH
  if (!is.numeric(object@stH)) {
    msg <- "'stH' should of type 'numeric'"
    errors <- c(errors, msg)
  }
  if (any(object@stH < 0)){
    msg <- "all 'stH' should be positiv or zero"
    errors <- c(errors, msg)
  }
  ## Lxh
  if (!is.numeric(object@Lxh)) {
    msg <- "'Lxh' should of type 'numeric'"
    errors <- c(errors, msg)
  }
  if (any(object@Lxh < 0)){
    msg <- "all 'Lxh' should be positiv or zero"
    errors <- c(errors, msg)
  }
  ## Hkz
  if (!is.integer(object@Hkz)) {
    msg <- "'Hkz' should of type 'integer'"
    errors <- c(errors, msg)
  }
  if (any(object@Hkz < 0) | any(object@Hkz > 2)){
    msg <- "all 'Hkz' should be 0, 1 or 2"
    errors <- c(errors, msg)
  }
  ## Skz
  if (!is.integer(object@Skz)) {
    msg <- "'Skz' should of type 'integer'"
    errors <- c(errors, msg)
  }
  if (any(object@Skz < 0) | any(object@Skz > 5)){
    msg <- "all 'Skz' should be inside 0:5"
    errors <- c(errors, msg)
  }
  ## Hsh
  if (!is.numeric(object@Hsh)) {
    msg <- "'Hsh' should of type 'numeric'"
    errors <- c(errors, msg)
  }
  if (any(object@Hsh < 0)){
    msg <- "all 'Hsh' should be positiv or zero"
    errors <- c(errors, msg)
  }
  ## Zsh
  if (!is.numeric(object@Zsh)) {
    msg <- "'Zsh' should of type 'numeric'"
    errors <- c(errors, msg)
  }
  if (any(object@Zsh < 0)){
    msg <- "all 'Zsh' should be positiv or zero"
    errors <- c(errors, msg)
  }
  ## Lsh
  if (!is.numeric(object@Lsh)) {
    msg <- "'Lsh' should of type 'numeric'"
    errors <- c(errors, msg)
  }
  if (any(object@Lsh < 0)){
    msg <- "all 'Lsh' should be positiv or zero"
    errors <- c(errors, msg)
  }
  ## Zab
  if (!is.numeric(object@Zab)) {
    msg <- "'Zab' should of type 'numeric'"
    errors <- c(errors, msg)
  }
  if (any(object@Zab < 0)){
    msg <- "all 'Zab' should be positiv or zero"
    errors <- c(errors, msg)
  }
  ## Lab
  if (!is.numeric(object@Lab)) {
    msg <- "'Lab' should of type 'numeric'"
    errors <- c(errors, msg)
  }
  if (any(object@Lab < 0)){
    msg <- "all 'Lab' should be positiv or zero"
    errors <- c(errors, msg)
  }
  ## Az
  if (!is.numeric(object@Az)) {
    msg <- "'Az' should of type 'numeric'"
    errors <- c(errors, msg)
  }
  if (any(object@Az < 0)){
    msg <- "all 'Az' should be positiv or zero"
    errors <- c(errors, msg)
  }
  ## LIh
  if (!is.numeric(object@LIh)) {
    msg <- "'LIh' should of type 'numeric'"
    errors <- c(errors, msg)
  }
  if (any(object@LIh < 0)){
    msg <- "all 'LIh' should be positiv or zero"
    errors <- c(errors, msg)
  }
  ## trL
  if (!is.numeric(object@trL)) {
    msg <- "'trL' should of type 'numeric'"
    errors <- c(errors, msg)
  }
  if (any(object@trL < 0)){
    msg <- "all 'trL' should be positiv or zero"
    errors <- c(errors, msg)
  }
  ## fixN
  if (!is.integer(object@fixN)) {
    msg <- "'fixN' should of type 'integer'"
    errors <- c(errors, msg)
  }
  if (any(object@fixN < 0)){
    msg <- "all 'fixN' should be positiv or zero"
    errors <- c(errors, msg)
  }
  ## fixL
  if (!is.numeric(object@fixL)) {
    msg <- "'fixL' should of type 'numeric'"
    errors <- c(errors, msg)
  }
  if (any(object@fixL < 0)){
    msg <- "all 'fixL' should be positiv or zero"
    errors <- c(errors, msg)
  }
  ## fixZ
  if (!is.numeric(object@fixZ)) {
    msg <- "'fixZ' should of type 'numeric'"
    errors <- c(errors, msg)
  }
  if (any(object@fixZ < 0)){
    msg <- "all 'fixZ' should be positiv or zero"
    errors <- c(errors, msg)
  }
  ## fixA
  if (!is.numeric(object@fixA)) {
    msg <- "'fixA' should of type 'numeric'"
    errors <- c(errors, msg)
  }
  if (any(object@fixA < 0)){
    msg <- "all 'fixA' should be positiv or zero"
    errors <- c(errors, msg)
  }
  ## fixR
  if (!is.numeric(object@fixR)) {
    msg <- "'Lxh' should of type 'numeric'"
    errors <- c(errors, msg)
  }
  if (any(object@fixR < 0) | any(object@fixR > 1)){
    msg <- "all 'fixR' should be inside [0,1]"
    errors <- c(errors, msg)
  }
  ## check parameter length
  lunms <- length(unique(sapply(slotNames(object), function(nms){
    length(slot(object, nms))
  }))) # length of unique slot length
  if(lunms != 1){
    msg <- "length of all slots should be equal."
    errors <- c(errors, msg)
  }

  if (length(errors) == 0) TRUE else errors
}

## set class ####
#' @title An S4 class to represent the parameters for tree assorting.
#' @description This class represents one or multiple parameter sets holding the
#' necessary information to specify the assortment process.
#' @slot stH stump height
#' @slot Lxh length of unusable wood at stem foot, see details
#' @slot Hkz height indicator, see details
#' @slot Skz stem indicator, see details
#' @slot Hsh height of stem wood, see details
#' @slot Zsh cutting diameter of stem wood, see details
#' @slot Lsh length of stem wood, see details
#' @slot Zab cutting diameter of upper trunk, see details
#' @slot Lab length of upper trunk, see details
#' @slot Az minimal cutting diameter, defaults to 7cm, see details
#' @slot LIh length of industrial wood, see details
#' @slot trL maximum transport length
#' @slot fixN number of fixed length assortments, see details
#' @slot fixL length of fixed length assortments, see details
#' @slot fixZ cutting diameter of fixed length assortments, see details
#' @slot fixA absolute add-on for good measure of fixed length assortments,
#' see details
#' @slot fixR relative add-on for good measure of fixed length assortments,
#' see details
#' @details The assortment process is defined by several parameters. These
#' follow the specification of its ancestor BDAT , but
#' are extended to allow for fix length assortments at the tree top
#' (industrial wood / pulp wood) and relaxes transport length and stump height.
#' \itemize{
#'   \item stH: stump height, defaults to 0, i.e. 1\% of tree height
#'   \item Lxh: length of unusable wood at stem foot [m], defaults to 0 (X-Holz)
#'   \item Hkz: indicator for tree top, 0 - normal, 1 - Wipfelbruch,
#'   2 - Gipfelbruch
#'   \itemize{
#'     \item 0 => H=H (default)
#'     \item 1 => H=H+2
#'     \item 2 => DBH < 30 => H=DBH; dbh > 30 => H = 30 + (DBH-30) * 0.3
#'   }
#'   \item Skz: indicator for stem type, defaults to 0
#'   \itemize{
#'     \item 0 => conifer trees => no assortment restriction;
#'     deciduous trees => no assortments
#'     \item 1 => monopodial deciduous trees => Hsh = 0.7*H
#'     \item 2 => branching between dbh and 7m => Hsh = 5m
#'     \item 3 => crown base < 3m => Hsh=0.1
#'     \item 4 => dead or broken stem => Az = H*0.7
#'     \item 5 => dead tree => non-usable wood
#'   }
#'   \item Hsh: usable stem height, defaults to 0, i.e. 0.7*H
#'   \item Zsh: minimum cutting diameter under bark for stem wood [cm], defaults
#'   to 0, using parameter \code{Az} if estimated length < maximum length
#'   (i.e. 20m)
#'   \item Lsh: length of stem wood, defaults to 0, i.e. length unrestricted
#'   \item Zab: minimum cutting diameter under bark for top segment [cm],
#'   defaults to 0, i.e. 14cm under bark
#'   \item Lab: length of top segment, defaults to 0, i.e. length unrestricted
#'   \item Az: minimum cutting diameter over bark [cm], defaults to 0,
#'   using an exponential function given DBH to estimate Az
#'   \item LIh: length of industrial wood [m], defaults to 0, i.e. length
#'   unrestricted
#'   \item trL: maximum transport length of assortments, defaults to 0, i.e. 19m
#'   \item fixN: number of fixed length assortments at stem foot, defaults to 0
#'   (no fixed length assortments, irrespective of other fix* parameters)
#'   \item fixZ: mininum diameter under bark for fixed length assortment at
#'   stem foot, defaults to 0
#'   \item fixL: length of fixed length assortment at stem foot, defaults to 0
#'   \item fixA: fixed length assortement add-on in [cm], defaults to 0
#'   \item fixR: fixed length assortement add-on in [\%], defaults to 0
#' }
#' @examples
#' parSort()
#' parSort(Lxh=1)
#' parSort(n=2)

setClass("parSort",
         representation(stH = "numeric",
                        Lxh  = "numeric",
                        Hkz  = "integer",
                        Skz  = "integer",
                        Hsh  = "numeric",
                        Zsh  = "numeric",
                        Lsh  = "numeric",
                        Zab  = "numeric",
                        Lab  = "numeric",
                        Az = "numeric",
                        LIh  = "numeric",
                        trL = "numeric",
                        fixN = "integer",
                        fixL = "numeric",
                        fixZ = "numeric",
                        fixA = "numeric",
                        fixR = "numeric"),
         prototype(stH=0, Lxh=0, Hkz=0L, Skz=0L, Hsh=0, Zsh=0, Lsh=0, Zab=14,
                   Lab=0, Az=0, LIh=0, trL=0, fixN=0L, fixL=0, fixZ=0, fixA=0,
                   fixR=0),
         validity = check_parSort)


## define subsetting ####
#' @title subsetting an object of class 'parSort'
#' @description using indices i and j to subset
## #' @name [-method
#' @describeIn parSort subsetting for class 'parSort'
#' @aliases [,parSort-method
#' @aliases [,parSort,ANY,ANY,ANY-method
#' @keywords methods
#' @docType methods
#' @param x object from which to extract
#' @param i index i
#' @param j index j
#' @param ... not currently used
#' @param drop drop dimensions, defaults to FALSE
#' @return a part of the original object
#' @export
setMethod("[", signature(x = "parSort"),
          function(x, i, j, ..., drop=FALSE){
            if(!is(x, "parSort"))
              stop("'x' needs to be of class 'parSort'")
            validObject(x)
            x@stH = x@stH[i]
            x@Lxh = x@Lxh[i]
            x@Hkz = x@Hkz[i]
            x@Skz = x@Skz[i]
            x@Hsh = x@Hsh[i]
            x@Zsh = x@Zsh[i]
            x@Lsh = x@Lsh[i]
            x@Zab = x@Zab[i]
            x@Lab = x@Lab[i]
            x@Az = x@Az[i]
            x@LIh = x@LIh[i]
            x@trL = x@trL[i]
            x@fixN = x@fixN[i]
            x@fixL = x@fixL[i]
            x@fixZ = x@fixZ[i]
            x@fixA = x@fixA[i]
            x@fixR = x@fixR[i]
            return(x)
          })

## constructor ####
#' @title constructor for class parSort
#' @description function to call \code{new()} on class \code{parSort}
#' @param n the number of parameter sets to generate, defaults to 1
#' @param stH stump height
#' @param Lxh length of unusable wood at stem foot, see details
#' @param Hkz height indicator, see details
#' @param Skz stem indicator, see details
#' @param Hsh height of stem wood, see details
#' @param Zsh cutting diameter of stem wood, see details
#' @param Lsh length of stem wood, see details
#' @param Zab cutting diameter of upper trunk, see details
#' @param Lab length of upper trunk, see details
#' @param Az minimal cutting diameter, defaults to 7cm, see details
#' @param LIh length of industrial wood, see details
#' @param trL maximum transport length
#' @param fixN number of fixed length assortments, see details
#' @param fixL length of fixed length assortments, see details
#' @param fixZ cutting diameter of fixed length assortments, see details
#' @param fixA absolute add-on for good measure of fixed length assortments,
#' given in cm; see details
#' @param fixR relative add-on for good measure of fixed length assortments,
#' given in percentage, i.e. 1\% = 1; see details
#' @param ... currently unused
#' @details if n is not given (or one) and any of the other parameter is given
#' with length greater than one, n is reset to the maximum length of all
#' parameters; care should be taken when using n and individual parameter
#' setting for several trees.
#' @return an object of class \code{parSort}, i.e. a list, each element of
#' length \code{n} or maximum of length of defined parameters
#' @export
parSort <- function(n=1, stH=0, Lxh=0, Hkz=0L, Skz=0L, Hsh=0, Zsh=0, Lsh=0,
                    Zab=14, Lab=0, Az=0, LIh=0, trL=0, fixN=0L, fixL=0, fixZ=0,
                    fixA=0, fixR=0, ...){
  ## eventually n is not given but one parameter has length bigger 1.
  ml <- max(c(length(stH), length(Lxh), length(Hkz), length(Skz), length(Hsh),
              length(Zsh), length(Lsh), length(Zab), length(Lab), length(Az),
              length(LIh), length(trL), length(fixN), length(fixL),
              length(fixZ), length(fixA), length(fixR)))
  if(n==1 & any(ml>1)) n <- ml

  new("parSort",
      stH = rep(stH, ifelse(length(stH) > 1, 1, n)),
      Lxh = rep(Lxh, ifelse(length(Lxh) > 1, 1, n)),
      Hkz = as.integer(rep(Hkz, ifelse(length(Hkz) > 1, 1, n))),
      Skz = as.integer(rep(Skz, ifelse(length(Skz) > 1, 1, n))),
      Hsh = rep(Hsh, ifelse(length(Hsh) > 1, 1, n)),
      Zsh = rep(Zsh, ifelse(length(Zsh) > 1, 1, n)),
      Lsh = rep(Lsh, ifelse(length(Lsh) > 1, 1, n)),
      Zab = rep(Zab, ifelse(length(Zab) > 1, 1, n)),
      Lab = rep(Lab, ifelse(length(Lab) > 1, 1, n)),
      Az = rep(Az, ifelse(length(Az) > 1, 1, n)),
      LIh = rep(LIh, ifelse(length(LIh) > 1, 1, n)),
      trL = rep(trL, ifelse(length(trL) > 1, 1, n)),
      fixN = as.integer(rep(fixN, ifelse(length(fixN) > 1, 1, n))),
      fixL = rep(fixL, ifelse(length(fixL) > 1, 1, n)),
      fixZ = rep(fixZ, ifelse(length(fixZ) > 1, 1, n)),
      fixA = rep(fixA, ifelse(length(fixA) > 1, 1, n)),
      fixR = rep(fixR, ifelse(length(fixR) > 1, 1, n)),
      ...)
}
