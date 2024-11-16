#' @title parameter of bark function
#' @description extract parameter of bark functions according to Altherr et al.
#' 1974 - 1979
#' @param ba tree species code; returned by \code{\link{BaMap}}
#' @param fn function number; see details
#' @param par parameter; see details
#' @details Function extracts the parameter according to tree species, function
#' type and parameter number. There are three parameters in each of four
#' functions. The first one refers to butt log (dt. Erdstamm), the second to
#' middle log (dt. Mittelstammstück), the third to the top log
#' (dt. Gipfelstammstück) and the fourth to the complete stem (dt. Gesamtstamm).
RiPar <- function(ba=NULL, fn=NULL, par=NULL){

  if(!is.null(ba)){
    if(ba < 1 | ba > 28)
      stop("'ba' must be integer inside interval [1,28]!")
  } else {
    ba <- 1:28
  }
  if(!is.null(fn)){
    if(fn < 1 | fn > 4)
      stop("'fn' must be integer inside interval [1,4]!")
  } else {
    fn <- 1:4
  }
  if(!is.null(par)){
    if(par < 1 | par > 3)
      stop("'par' must be integer inside interval [1,3]!")
  } else {
    par <- 1:3
  }

  return(dbtp[ba, fn, par])
}
