# TapeS 0.13.1-9000

* in tprVolume now checking lower end 'A' of segment against upper end 'B' of
  segment: A should be lower or equal B. This is especially important if one of
  both is given as diameter and the other as a height.
* bug fix in uncertainty estimation of biomass for deciduous trees and requested
  component 'ndl'.

# TapeS 0.13.1

* extended internal function which calculated position and length of measurement
  position during volume estimation: now parameter AB in tprVolume can both be
  zero with no error and returning reasonable information: 0.
  
# TapeS 0.13.0

* implemented variance estimates for biomass estimates, based on the NSUR 
  component biomass functions. See details in tprBiomass-helpfile for how this 
  was done.
* added function nsur2() implementing the 'simple' biomass functions of 
  Vonderach et al. 2018 (internal).
* corrected function nsur(), now using the adjusted tree height (tree height minus
  stump height) for prediction as detailed in Vonderach & Kändler (2021)
* bugfix in internal nsur()-function: one parameter in ash was wrong (led to
  negative estimates in small trees).
* extended tprBiomass: now optionally scale NSUR- to NFI-estimate
* moved internal data to sysdata instead of data/, hence, package runs also if
  not attached (e.g. calls via TapeS::tprVolume don't break)
* added citation file
* corrected minor bugs, some internal improvements

# TapeS 0.12.1

* corrected bug in example of nsur() function

# TapeS 0.12.0

* function tprVolume now with MSE and interval information
* function tprDiameter now returns estimated MSE for the mean and prediction
  when interval='MSE'
* added 'sHt' into tprObject
* added function FormTariff() to fix taper form according to German NFI
  observations instead of giving a second measured diameter
* added 'inv'-parameter into tprTrees-constructor
* corrected component biomass function, now total mass includes needle mass
* improved plotting function
* improved documentation


# TapeS 0.11.2

* updated the relative taper models for beech, oak and red oak, adding a
  lower diameter to the relative data based on BDAT to stabilise taper curves.

# TapeS 0.11.1

* small adjustements

# TapeS 0.11.0

* added new TapeR-models for beech, oak and red oak based on pre-smoothed
  data used during BDAT model building. These models now are the default models
  for these species; the original TapeR-models based one sectional data are now
  pushed to species code 37, 38 and 39, respectively. Updated all functions to
  handle the extended data set.
* corrected bug in function HDxoR()
* technical adjustments

# TapeS 0.10.0

* added height tariff based on German NFI data to estimated tree height
  given diameter in breast height (DBH)

# TapeS 0.9.0

* corrected bug in bdat_as_tprtrees(), when passing more than one tree from
  rBDAT::buildTree (class 'datBDAT')
* added 37th tree species based on an additional taper curve model, which is
  build on pre-smoothed beech data used during development of BDAT.
* corrected bug in plot function
* added function and indicator for residual variance to be set to zero to
  force taper curve through measurements; even more options possible
* added support diameter D001 (in 1\% of height) at the stem base to avoid
  non-monotoneous taper curves for such trees (mainly very small trees of
  certain species). The support diameter is added based on the population
  average relation between D001 and D05.

# TapeS 0.8.x

* corrected bug in tprDiameter, when bark=FALSE and interval != 'none'
* corrected biomass estimation: tprBiomass now return aboveground biomass
  *without* component needles as the underlying functions do not include 
  needles.  

# TapeS 0.8.5

* more flexibility in defining a tprTrees object
* corrected bug, when variable HdhG was not always available in tprAssortment
* corrected bug, when indexing was not correct deciduous trees in tprAssortment

# TapeS 0.8.4

* simplified the return value of tprDiameterCpp (as in tprDiameter)
* added examples

# TapeS 0.8.3

* allow for small diameter increase in taper curve when checking monotonicity
  (increase <= 1% of diameter)
* constraint in bark function: never return negative double bark thickness 
  (c.f. e.g. Japanese larch) and diameter after bark deduction might not be 
  smaller than zero

# TapeS 0.8.2

* update examples in tprDiameter and tprDiameterCpp

# TapeS 0.8.1

* updated helpfile of class parSort

# TapeS 0.8.0

* plot method for class 'tprTrees'
* call of 'tprHeight' in 'tprAssortment' corrected
* in subsetting: keep attribute about R0-status in object

# TapeS 0.7.1

* updated helpfile for tprVolume

# TapeS 0.7.0

* added confidence and prediction intervals into tprVolume
* added simTrees()
* some fixes in show and nfi_as_tprtrees coercion

# TapeS 0.6.0

* renamed 'trees'-class into 'tprTrees' to be more specific. 
* renamed slot 'H' into 'Ht' (and accessort function)
* added cartesian product to function tprDiameter/Cpp, tprHeight and tprBark
* some updates to BDAT-coercion and show-method
* some updates to vignette

# TapeS 0.5.1

* switched option TapeS_R0 from TRUE to FALSE (as default)

# TapeS 0.5.0

* added accessor functions to objects of class 'trees'
* added 'show' function

# TapeS 0.4.0

* re-implemented component biomass estimation to get rid of 
  dplyr/tidyr-dependency. Now component biomass is estimated using newly
  developed NSUR-biomass functions. Implemented using Rcpp.
* added coercion function to coerce nfi, sectional and bdat data
  
# TapeS 0.3.0

* removed dependency on (unpublished) rBDATPRO-package, implementing biomass 
  function via RCpp (instead call to rDBATPRO)
* added monotonicity check of taper curve right into class-constructor

# TapeS 0.2.0

* included fitted TapeR models (data: SKPar) and updated corresponding data 
  help file. 
  Added this News file.

# TapeS 0.1.0

* basic package development



