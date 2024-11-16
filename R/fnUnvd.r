#' @title percentage of unusable coarse wood
#' @description function extracts the percentage of unusable coarse wood
#' according to species (beech, oak), diameter class and cutting diameter
#' @param ba tree species index; see details
#' @param dm diameter class; see details
#' @param cd cutting diameter; see details
#' @details Function extracts the percentage of unusable coarse wood according
#' to three parameters: (i) tree species, which is 1 for using beech models and
#' 2 for using the oak model; (ii) the 2cm-diameter class (from 8 and 60cm) and
#' (iii) the cutting diameter ranging from 8 to 40cm.
#' @references Kublin and Scharnagl (1988): Verfahrens- und Programmbeschreibung
#' zum BWI-Unterprogramm BDAT. FVA-BW 1988. ISSN: 0178-3165.

fnUnvd <- function(ba=NULL, dm=NULL, cd=NULL){

  if(!is.null(ba)){
    if(ba < 1 | ba > 2)
      stop("'ba' must be integer inside interval [1,2]!")
  } else {
    ba <- 1:2
  }
  if(!is.null(dm)){
    if(dm < 1 | dm > 27)
      stop("'dm' must be integer inside interval [1,27]!")
  } else {
    fn <- 1:27
  }
  if(!is.null(cd)){
    if(cd < 1 | cd > 33)
      stop("'cd' must be integer inside interval [1,33]!")
  } else {
    par <- 1:33
  }
  #
  return(unvd[ba, dm, cd])
}

