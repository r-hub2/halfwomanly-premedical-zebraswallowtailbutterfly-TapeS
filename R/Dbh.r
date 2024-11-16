#' @title Extract pre-defined diameter
#' @description Function extracts pre-defined diameters, e.g. dbh (in 1.3m) or
#' D03 (in 30\% of tree height) for a \code{\link{tprTrees}}-object
#' @param obj a object of class 'tprTrees'
#' @details a wrapper around \code{\link{tprDiameter}} to calculate specifically
#' defined diameters like diameter in breast height (dbh), diameter in 7m above
#' ground or in 5\% and 30\% of tree height.
#' @return diameter(s) in predefined heights
#' @name Dbh
#' @examples
#' t <- tprTrees()
#' Dbh(t) # diameter in breast height (i.e. 1.3m)
#' Bhd(t) # same, german named function name
#' D13(t) # same, height related function name
#' D005(t) # diameter in 5% of tree height
#' D7(t) # diameter in height of 7m
#' D03(t) # diameter in 30% of tree height
NULL
#' @describeIn Dbh wrapper to calculate diameter in breast height
#' @export
Dbh <- function(obj){
  D <- tprDiameterCpp(obj, Hx=1.3, bark=TRUE)
  return(D)
}
#' @describeIn Dbh German alias for function Dbh
#' @export
Bhd <- Dbh
#' @describeIn Dbh Height specific alias for function Dbh
#' @export
D13 <- Dbh

#' @describeIn Dbh Function to calculate diameter over bark in 7m above ground
#' @export
D7 <- function(obj){
  D <- tprDiameterCpp(obj, Hx=7, bark=TRUE)
  return(D)
}

#' @describeIn Dbh Function to calculate diameter over bark in 30\% of tree height
#' @export
D03 <- function(obj){
  D <- tprDiameterCpp(obj, Hx=0.3*Ht(obj), bark=TRUE, cp=FALSE)
  return(D)
}

#' @describeIn Dbh Function to calculate diameter over bark in 5\% of tree height
#' @export
D005 <- function(obj){
  D <- tprDiameterCpp(obj, Hx=0.05*Ht(obj), bark=TRUE, cp=FALSE)
  return(D)
}
