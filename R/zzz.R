#'@import tidyverse
#'@import RTMB
#'@import memoise

.onLoad <- function(libname, pkgname) {
	options(show.error.locations = TRUE)
	options(keep.source = TRUE)
	options(keep.source.pkgs = TRUE)

	options(warnPartialMatchDollar = TRUE)
	options(warnPartialMatchAttr = TRUE)
	options(warnPartialMatchArgs = FALSE)

	#cat("Here 123!\n")
	#print("We are here")
	
	#library(memoise)
	environment(memoise_seq_bastl) <<- asNamespace("memoise")
	#environment(memoise_seq_bastl) <<- new.env(parent = asNamespace("memoise"))
}
