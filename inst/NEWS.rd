\name{NEWS}
\title{NEWS}

\section{Changes in version 0.12.1}{
\itemize{
\item corrected bug in example of nsur() function
\item added CITATION file
}
}

\section{Changes in version 0.12.0}{
\itemize{
\item function tprVolume now with MSE and interval information
\item function tprDiameter now returns estimated MSE for the mean and prediction
  when interval='MSE'
\item added 'sHt' into tprObject
\item added function FormTariff() to fix taper form according to German NFI
  observations instead of giving a second measured diameter
\item added 'inv'-parameter into tprTrees-constructor
\item corrected component biomass function, now total mass includes needle mass
\item improved plotting function
\item improved documentation
}
}

\section{Changes in version 0.11.2}{
\itemize{
\item updated the relative taper models for beech, oak and red oak, adding a
  lower diameter to the relative data based on BDAT to stabilise taper curves.
}
}

\section{Changes in version 0.11.1}{
\itemize{
\item small adjustements
}
}

\section{Changes in version 0.11.0}{
\itemize{
\item added new TapeR-models for beech, oak and red oak based on pre-smoothed
  data used during BDAT model building. These models now are the default models
  for these species; the original TapeR-models based one sectional data are now
  pushed to species code 37, 38 and 39, respectively. Updated all functions to
  handle the extended data set.
\item corrected bug in function HDxoR()
\item technical adjustments
}
}

\section{Changes in version 0.10.0}{
\itemize{
\item added height tariff based on German NFI data to estimated tree height
  given diameter in breast height (DBH)
}
}

\section{Changes in version 0.9.0}{
\itemize{
\item corrected bug in bdat_as_tprtrees(), when passing more than one tree from
  rBDAT::buildTree (class 'datBDAT')
\item added 37th tree species based on an additional taper curve model, which is
  build on pre-smoothed beech data used during development of BDAT.
\item corrected bug in plot function
\item added function and indicator for residual variance to be set to zero to
  force taper curve through measurements; even more options possible
\item added support diameter D001 (in 1\% of height) at the stem base to avoid
  non-monotoneous taper curves for such trees (mainly very small trees of
  certain species). The support diameter is added based on the population
  average relation between D001 and D05.
}
}

\section{Changes in version 0.8.6}{
\itemize{
\item corrected bug in tprDiameter, when bark=FALSE and interval != 'none'
\item corrected biomass estimation: tprBiomass now returns aboveground biomass
 *without* component needles as the underlying functions do not include needles.
}
}

\section{Changes in version 0.8.5}{
\itemize{
\item more flexibility in defining a tprTrees object
\item corrected bug, when variable HdhG was not always available in tprAssortment
\item corrected bug, when indexing was not correct deciduous trees in tprAssortment
}
}

\section{Changes in version 0.8.4}{
\itemize{
\item simplified the return value of tprDiameterCpp (as in tprDiameter)
\item added examples
}
}

\section{Changes in version 0.8.3}{
\itemize{
\item allow for small diameter increase in taper curve when checking monotonicity
(increase <= 1\% of diameter)
\item constraint in bark function: never return negative double bark thickness
(c.f. e.g. Japanese larch) and diameter after bark deduction might not be
smaller than zero
}
}

\section{Changes in version 0.8.2}{
\itemize{
\item update examples in tprDiameter and tprDiameterCpp
}
}

\section{Changes in version 0.8.1}{
\itemize{
\item updated helpfile of class parSort
}
}

\section{Changes in version 0.8.0}{
\itemize{
\item plot method for class 'tprTrees'
\item call of 'tprHeight' in 'tprAssortment' corrected
\item in subsetting: keep attribute about R0-status in object
}
}

\section{Changes in version 0.7.1}{
\itemize{
\item updated helpfile for tprVolume
}
}

\section{Changes in version 0.7.0}{
\itemize{
\item added confidence and prediction intervals into tprVolume
\item added simTrees()
\item some fixes in show and nfi_as_tprtrees coercion
}
}

\section{Changes in version 0.6.0}{
\itemize{
\item renamed 'trees'-class into 'tprTrees' to be more specific.
\item renamed slot 'H' into 'Ht' (and accessort function)
\item added cartesian product to function tprDiameter/Cpp, tprHeight and tprBark
\item some updates to BDAT-coercion and show-method
\item some updates to vignette
}
}

\section{Changes in version 0.5.1}{
\itemize{
\item switched option TapeS_R0 from TRUE to FALSE (as default)
}
}

\section{Changes in version 0.5.0}{
\itemize{
\item added accessor functions to objects of class 'trees'
\item added 'show' function
}
}

\section{Changes in version 0.4.0}{
\itemize{
\item re-implemented component biomass estimation to get rid of
dplyr/tidyr-dependency. Now component biomass is estimated using newly
developed NSUR-biomass functions. Implemented using Rcpp.
\item added coercion function to coerce nfi, sectional and bdat data
}
}

\section{Changes in version 0.3.0}{
\itemize{
\item removed dependency on (unpublished) rBDATPRO-package, implementing biomass
function via RCpp (instead call to rDBATPRO)
\item added monotonicity check of taper curve right into class-constructor
}
}

\section{Changes in version 0.2.0}{
\itemize{
\item included fitted TapeR models (data: SKPar) and updated corresponding data
help file.
Added this News file.
}
}

\section{Changes in version 0.1.0}{
\itemize{
\item basic package development
}
}

