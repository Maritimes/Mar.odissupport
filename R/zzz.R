.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste0("Version: ", utils::packageDescription('Mar.odissupport')$Version))
}
.onLoad <- function(libname, pkgname){
  options(stringsAsFactors = FALSE)
}
