.onLoad <- function(libname, pkgname){
  setTapeSoptions(Rfn = list(fn="sig2"), mono = TRUE)
}

#' @importFrom utils packageVersion
.onAttach <- function(libname, pkgname){
  setTapeSoptions(Rfn = list(fn="sig2"), mono = TRUE)
  optMessage <- paste0("Option 'TapeS_Rfn' set to '", options()$TapeS_Rfn, "'.",
                       " See '?setTapeSoptions' for details.")
  packageStartupMessage(paste0("This is TapeS ", packageVersion("TapeS"), "\n",
                               optMessage))
}

.onUnload <- function (libpath) {
  library.dynam.unload("TapeS", libpath)
}
