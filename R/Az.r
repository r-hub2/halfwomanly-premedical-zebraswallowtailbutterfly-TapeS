#' @title estimate minimum cutting diameter
#' @description estimate minimum cutting diameter based on diameter in breast
#' height based on the functions implemented in BDAT
#' @param sp Bdat species code [1;36], integer
#' @param dbh vector of diameter in breast height, numeric
#' @details the implemented BDAT function and parameters are used. Not all
#' BDAT-species possess their own parameters, hence most of them are matched to
#' one of the main tree species, especially in deciduous tree species (only
#' parameters for beech and oak are available).
#' @return vector of minimum cutting diameter [cm].
#' @examples
#' sp <- 1
#' dbh <- 30
#' Az(sp, dbh)
#' @export

Az <- function(sp, dbh){
  if(length(sp) > 1 & length(dbh) > 1 & length(sp) != length(dbh)){
    sp <- sp[1]
    warning("length of sp and dbh unequal: using first element of sp only!")
  }
  if(any(sp > 36) | any(sp < 1)){
    sp <- 1
    warning("species code reset to 1 (spruce)!")
  }
  if(!is.numeric(dbh) | any(dbh <= 0)) stop("dbh must be a positive numeric!")

  sp <- BaMap(sp, 3) # remapping species code by using specific function

  retAz <- exp(azp[sp, 1] +
               azp[sp, 2] * log(dbh) +
               azp[sp, 3] * (log(dbh)^2) +
               azp[sp, 4] * dbh)
  return(retAz)
}
