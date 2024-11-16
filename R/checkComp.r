#' generate and/or check validity of biomass function component names
#'
#' @param component vector of biomass component names, see details
#' @details If \code{component} is NULL, by default, component name
#' for total aboveground biomass is returned. If is \code{all},
#' then all available component names are returned.
#' \itemize{
#'  \item{stw: stump wood}
#'  \item{stb: stump bark}
#'  \item{sw: solid wood with diameter above 7cm over bark}
#'  \item{sb: bark of component sw}
#'  \item{fwb: fine wood incl. bark}
#'  \item{ndl: needles}
#'  \item{agb: total aboveground biomass}
#' }
#' @return a vector of component names
#' @examples
#' \dontrun{
#' TapeS:::checkComp()
#' TapeS:::checkComp("AGB")
#' TapeS:::checkComp("biomass")
#' }
#'
check_Comp <- function(component=NULL){
  if(is.null(component)){
    comp <- "agb"
  } else if(identical(tolower(component), "all")){
    comp <- c("stw", "stb", "sw", "sb", "fwb", "ndl", "agb")
  } else if(!all(tolower(component) %in% c("stw", "stb", "sw", "sb",
                                           "fwb", "ndl", "agb"))){
    stop("'component' wrong! use 'stw', 'stb', 'sw', 'sb', ",
         "'fwb', 'ndl' and/or 'agb'")
  } else {
    comp <- tolower(unique(component))
  }
  return(comp)
}
