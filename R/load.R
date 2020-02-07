.onLoad <- function(libname, pkgname) {
  # print("onLoad")
}


.onUnload <- function (libpath) {
  cleanup()
  # print("onUnload")
  library.dynam.unload("ParallelSampEn", libpath)
}
