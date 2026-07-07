# E.g. when making source tarballs, do it like this
# Rscript tools/build.R
# rather than bare R CMD build !

source("tools/write-build-info.R")
system2(file.path(R.home("bin"), "R"), c("CMD", "build", "."))
