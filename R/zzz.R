.First.lib <- function(lib, pkg) library.dynam("mscalib",pkg,lib)
.Last.lib <- function(libpath) library.dynam.unload("mscalib", libpath)

