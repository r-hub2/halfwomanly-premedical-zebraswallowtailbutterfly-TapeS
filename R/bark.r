#' @title Double Bark Thickness
#' @description Function returns double bark thickness according to
#' Altherr et al. 1974/75/76/78/79
#' @param Ba tree species according to BDAT, cf. \code{\link{tprSpeciesCode}}
#' @param Dm diameter for which double bark thickness is requested
#' @param relH relative height of Dm inside stem
#' @details Function re-implemented according to
#' Subroutine RINDE(Hhrel,Kw,Ri,Hsga,Zo), BDAT-fortran Code line 5691ff.
#' No Functions for (historic) Heilbronner Sortierung implemented.
#'
#' NB: to avoid negative double bark thickness, such values are constraint to
#' zero. Additionally, diameter after bark reduction might not be smaller than
#' zero, hence double bark thickness is reduce to \code{Dm}.
#' @return double bark thickness [cm]
#' @references
#' Altherr, E., P. Unfried, J. Hradetzky and V. Hradetzky (1974). Statistische
#' Rindenbeziehungen als Hilfsmittel zur Ausformung und Aufmessung unentrindeten
#' Stammholzes. Kiefer, Buche, Hainbuche, Esche und Roterle. Freiburg i. Br.,
#' Forstl. Versuchs- u. Forschungsanst. Baden-Württenberg.
#'
#' Altherr, E., P. Unfried, J. Hradetzky and V. Hradetzky (1975). Statistische
#' Rindenbeziehungen als Hilfsmittel zur Ausformung und Aufmessung unentrindeten
#' Stammholzes. Europäische Lärche, Japanische Lärche, Schwarzkiefer, Stieleiche,
#' Traubeneiche, Roteiche, Bergahorn und Linde. Freiburg i. Br., Forstl.
#' Versuchs- u. Forschungsanst. Baden-Württenberg.
#'
#' Altherr, E., P. Unfried, J. Hradetzky and V. Hradetzky (1976). Statistische
#' Rindenbeziehungen als Hilfsmittel zur Ausformung und Aufmessung unentrindeten
#' Stammholzes. Weymouthskiefer, Robinie, Bergulme, Birke, Marilandica-Pappel
#' und Robusta-Pappel. Freiburg i. Br., Forstl. Versuchs- u. Forschungsanst.
#' Baden-Württenberg.
#'
#' Altherr, E., P. Unfried, J. Hradetzky and V. Hradetzky (1978). Statistische
#' Rindenbeziehungen als Hilfsmittel zur Ausformung und Aufmessung unentrindeten
#' Stammholzes. Fichte, Tanne, Douglasie und Sitka-Fichte. Freiburg i. Br.,
#' Forstl. Versuchs- u. Forschungsanst. Baden-Württenberg.
#'
#' Altherr, E., P. Unfried, J. Hradetzky and V. Hradetzky (1979). Statistische
#' Rindenbeziehungen als Hilfsmittel zur Ausformung und Aufmessung unentrindeten
#' Stammholzes. Neupotz-Pappel, Regenerata-Pappel, Kirsche, Spitzahorn, Feldahorn,
#' Aspe, Weide, Flatterulme, Tulpenbaum u. Elsbeere. Freiburg i. Br., Forstl.
#' Versuchs- u. Forschungsanst. Baden-Württenberg.
#' @export
#' @examples
#' bark(1, 30, .1)
#' bark(11, 4, .1) # zero instead of -0.2497
bark <- function(Ba, Dm, relH){

  stopifnot(identical(length(Ba), length(Dm)))
  stopifnot(identical(length(Ba), length(relH)))
  i <- ifelse((1-relH) <= 0.4, 3, # function for top log
              ifelse((1-relH) <= 0.7, 2, # function for middle log
                     1)) # function for butt log

  BaBark <- BaMap(Ba, 2)
  par <- sapply(seq(along=i), function(a){
    RiPar(ba = BaBark[a], fn = i[a])
  })
  DBT <- sapply(seq(along=Dm), function(a){
    #double bark thickness in cm
    (par[1,a] + par[2,a]*Dm[a] + par[3,a]*Dm[a]*Dm[a])*.1
  })
  # avoid negative double bark thickness (cf. Japanese larch)
  DBT <- ifelse(DBT < 0, 0, DBT)
  # avoid negative diameter under bark after bark reduction
  # in this case bark needs to be equal to dm
  DBT <- ifelse((Dm - DBT) < 0, Dm, DBT)


  return(DBT)
}
