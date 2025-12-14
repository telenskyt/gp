#@xxixmportFrom dplyr select filter arrange
#'@import dplyr
#'@importFrom purrr map
#'@import RTMB
#'@import memoise

# toto ne, pouzil jsem usethis::use_pipe()! viz https://stackoverflow.com/a/52231630/684229
# @xxximportFrom magrittr %>%

.onLoad <- function(libname, pkgname) {
	options(show.error.locations = TRUE)
	options(keep.source = TRUE)
	options(keep.source.pkgs = TRUE)
	options(pillar.print_max = 120)
	
	options(warnPartialMatchDollar = TRUE)
	options(warnPartialMatchAttr = TRUE)
	options(warnPartialMatchArgs = FALSE)

	#cat("Here 123!\n")
	#print("We are here")
	
	#library(memoise)
	environment(memoise_seq_bastl) <<- asNamespace("memoise")
	#environment(memoise_seq_bastl) <<- new.env(parent = asNamespace("memoise"))
}
