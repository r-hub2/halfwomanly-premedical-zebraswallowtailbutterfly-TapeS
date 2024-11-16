#' @title Set and get options for the TapeS-package
#' @description Function to set and get options on how the TapeS-package works.
#' @param Rfn setting for residuals error matrix, defaults to \code{"sig2"}, see
#'   details.
#' @param mono logical, defaults to true. If calibrated taper curve is
#'   non-monotonic at stem base, an support diameter is added.
#' @details So far, only two options are implemented: \code{TapeS_Rfn} and
#'   \code{TapeS_mono}. Teh first defaults to "sig2" (i.e. 'sigma squared') and
#'   the second to "TRUE".
#'
#'   The TapeR-taper curves can be evaluated in basically two ways: (i) either
#'   as defined in the TapeR-package, i.e. the diameters and volumes are
#'   estimated using the estimated error structure and find an optimal taper
#'   curve given the measured diameters or (ii) by interpolating the measured
#'   diameters, i.e. forcing the estimated taper curve through those
#'   measurements by setting the residual error structure to zero. See Kublin et
#'   al. (2013), p.987 bottom left. Technically, forcing the taper curve through
#'   the measurements is achieved by setting the residual error matrix R to
#'   zero, that is \code{Rfn = list(fn="zero")}. Defaults to \code{Rfn =
#'   list(fn="sig2")}. Besides, one can defined other functions about
#'   assumptions about the errors at the measurement positions, see
#'   \code{\link[TapeR]{resVar}} for options.
#'
#'   NB: Caution is required in applying \code{Rfn=list(fn="zero")}, since
#'   forcing the taper curve through too many points might lead to singularities
#'   or implausible results!
#'
#'   The option 'mono=TRUE' assures that no taper curve is generated which shows
#'   lower diameter in lower heights, possibly adding a support diameter at 1\%
#'   of tree height.
#'
#' @references Kublin, E., Breidenbach, J., Kaendler, G. (2013) A flexible stem
#'   taper and volume prediction method based on mixed-effects B-spline
#'   regression, Eur J For Res, 132:983-997.
#' @return by defaults, sets \code{options()$TapeS_Rfn} to "sig2"
#' @export
#' @examples
#' ## reset option TapeS_Rfn to "sig2", i.e. model based errors by
#' setTapeSoptions(Rfn = list(fn="sig2"))
#' ## or to force the taper curve through the  measurements, set
#' options("TapeS_Rfn" = list(fn="zero"))

setTapeSoptions <- function(Rfn = list(fn="sig2"), mono = TRUE){
  options(TapeS_Rfn = Rfn, TapeS_mono = mono)
}

#' @param name name of options to be returned
#' @rdname setTapeSoptions
#' @export
#' @examples
#' ## see the actual state of options by
#' options()[grep("^TapeS_", names(options()))]
#' ## or easier
#' getTapeSoptions()

getTapeSoptions <- function(name=NULL){
  r <- options()[grep("TapeS", names(options()))]
  if(is.numeric(name) | is.integer(name)){
    name <- as.integer(name)
    name <- name[name>0 & name < 2]
    name <- c("TapeS_mono", "TapeS_Rfn")[name]
  }
  if(any(!(name %in% c("TapeS_mono", "TapeS_Rfn")))){
    warning("name must only be 'TapeS_mono' or 'TapeS_Rfn'.")
  }
  name <- name[name %in% c("TapeS_mono", "TapeS_Rfn")]
  if(is.null(name)){
    return(r)
  } else {
    return(r[name])
  }

}
