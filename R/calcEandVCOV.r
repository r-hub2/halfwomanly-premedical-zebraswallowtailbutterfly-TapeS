#' @title calculate VCOV-Matrix for volume segments
#' @description calculate variance-covariance matrix for volume segments from
#' the estimated diameter and uncertainty information from TapeR-model
#' @param estD vector of estimated diameter, numeric
#' @param kovD variance-covariance-matrix of the estimated diameter, numeric
#' @param estL vector of segment length, numeric
#' @details Calculations according to rules for products and sums of variances
#' @return variance-covariance matrix of the segment volume
#' @importFrom stats qt cov2cor

calcVCOVsekVol <- function(estD, kovD, estL){
  Cm <- pi/4 * 1e-4 * estL # constant factor per segment
  varD <- diag(kovD)
  sdD <- sqrt(varD)
  corD <- cov2cor(kovD)

  # variance-covariance of segment volume m³
  # es gilt: Var(a*X) = a^2*Var(X)
  # daher: Var(Cm*D^2)=Cm^2*Var(X^2)=Cm^2 * (3*varD^2 + 4*varD * D^2)
  # vereinfacht aus var(x1*x1)= var(x1)*var(x1) + 2*cov(x1, x1)^2 + 4*mean(x1)*mean(x1)*cov(x1, x1)
  varSek <- Cm^2 * (3*diag(kovD)^2 + 4*diag(kovD)*estD^2)
  ## es wird die kovarianz für d1^2*l1, d2^2*l2, ..... benötigt
  ## cov(cm*x^2, cm*y^2) = cm^2 * (4*mean(x)*mean(y)*cor(x,y)*sd(x)*sd(y) + 2*cor(x,y)^2*var(x)*var(y))
  ## cf. last formula in appendix of Kublin et al (2013, p. 996)
  kovSek <- kovD # just copy for dimension
  diag(kovSek) <- Cm^2 * (3*diag(kovD)^2 + 4*diag(kovD)*estD^2)
  for(i in seq(along=estD)){
    for(j in seq(along=estD)){
      if(i>j){
        #i <- 2; j <- 1
        kovSek[i, j] <- Cm[i]*Cm[j] * (4*estD[i]*estD[j]*corD[i,j]*sdD[i]*sdD[j] + 2*corD[i,j]^2*varD[i]*varD[j])
        kovSek[j, i] <- kovSek[i, j]
      }
    }
  }
  return(kovSek)
}
