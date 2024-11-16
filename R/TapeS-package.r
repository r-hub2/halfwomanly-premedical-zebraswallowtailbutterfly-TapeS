#' TapeS: New taper functions for most common tree species in Germany and
#'   evaluating functions.
#'
#' Implementing new taper functions for German NFI and making available
#' associated function to evaluate diameter, height of diameter, bark thickness,
#' biomass and volume of these taper functions. Also providing wrappers to
#' calculate specific diameters like dbh, d005, d03 and d7 as well as specific
#' volume measures like Vfm, Efm, VolR, VolE, Vol_FAO, and physical stem volume.
#'
#' The package contains datasets which hold the taper curve models and further
#' parameters used inside the provided functions
#' (see \code{data(package="TapeS")}).
#'
#' @name TapeS-package
#' @aliases TapeS-package
#' @docType package
#' @importFrom utils globalVariables
#' @keywords package
utils::globalVariables(c("SKPar", "azp", "dbtp", "unvd", "volfaodlt7",
                         "crownExpansion"))
