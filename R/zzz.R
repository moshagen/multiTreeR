.onLoad <- function(libname, pkgname) {
  #  rJava::.jpackage(pkgname, lib.loc = libname)
  .jpackage(pkgname, lib.loc = libname)
}
