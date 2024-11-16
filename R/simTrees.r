#' @title simulating objects of class \code{tprTrees}
#' @description Function to simulate an object of class \code{tprTrees}
#' @param par list of lists, one for each species
#' @details Function simulates trees based on given distributions and petterson
#' height function. Dbh can be simulated using normal ('norm'), weibull or
#' gamma distribution. Others might be added.
#'
#' The \code{par}-list of each species needs the following named entries:
#' \code{spp} - species code, \code{n} - number of trees, \code{ddist} -
#' distribution of dbh, \code{dpar} - list of parameter of the distribution,
#' i.e. \code{mu} and \code{sd} for normal distribution and \code{shape} and
#' \code{scale} for weibull and gamma distribution. The latter both might use
#' \code{lag} to offset the estimated diameter by this amount.
#' @return an object of class \code{tprTrees}
#' @seealso \code{\link{petterson}} for the implemented height function and
#' \code{\link{dnorm}}, \code{\link{dweibull}} and \code{\link{dgamma}} for the
#' diameter distributions.
#' @importFrom stats rnorm rweibull rgamma
#' @export
#' @examples
#' par <- list(list(spp=1, n=10, ddist="norm", dpar=list(mu=30, sd=4)),
#'             list(spp=3, n=5, ddist="norm", dpar=list(mu=40, sd=2)))
#' simTrees(par)
#'
simTrees <- function(par=NULL){
  stopifnot(is.list(par) | is.null(par))
  if(is.null(par)){
    par <- list(list(spp=1, n=10, ddist="norm", dpar=list(mu=30, sd=4)))
  }
  ll <- list()
  for(i in seq(along=par)){
    pari <- par[[i]]
    n <- pari$n
    id <- seq(1, n)
    lag <- ifelse(!is.null(pari$dpar$lag), pari$dpar$lag, 0)
    if(pari$ddist == "norm"){
      d1 <- rnorm(n, mean = pari$dpar$mu, sd = pari$dpar$sd)
    } else if(pari$ddist == "weibull"){
      d1 <- rweibull(n, shape = pari$dpar$sh, scale = pari$dpar$sc) + lag
    } else if(pari$ddist == "gamma"){
      d1 <- rgamma(n, shape = pari$dpar$sh, scale = pari$dpar$sc) + lag
    } else {
      stop("'", pari$ddist, "' not implemented!", sep = "")
    }

    q03 <- 0.8
    Ht <- sapply(id, function(a) petterson(pari$spp, d1[a]))
    ll[[i]] <- data.frame(id = id,
                          BaTpr = pari$spp,
                          Bhd = d1 * 10,
                          H1 = 1.3 * 10,
                          D03 = d1 * q03 * 10,
                          H2 = 0.3 * Ht * 10,
                          Hoehe = Ht * 10)
  }
  res <- do.call(rbind, ll)
  res <- nfi_as_tprtrees(res)# nfi_as_tprtrees requires [mm] and [dm]
  return(res)
}

