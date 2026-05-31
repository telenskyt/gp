


Fisher_dispatcher <- function (gp, f, data = gp$data, hyperpar, fisher = NULL)
{
	stopifnot(!is.null(fisher))
	if (!is(gp[["lik"]], 'likTempPhased')) # !!!! ??? maybe not for the sampling???
		stop("For the Laplace-Fisher method, you need likTempPhased likelihood tempate, see ?gp.")
	if (!is.list(fisher))
		stop("parameter 'fisher' has to be a list")
	if (!is.null(fisher[["sampling"]]) && fisher$sampling == TRUE) {
		return(do.call("Fisher_sampling", list(gp, f, data, hyperpar, fisher)))
	} else {
		if (!is.character(gp$lik$likType))
			stop("For the Laplace-Fisher method, lik part of the likTempPhased template must be a character")
		if (gp$lik$likType == "bern")
			return(do.call("Fisher_bern", list(gp, f, data, hyperpar, fisher)))
		else
			stop("gp(lik$likType = '", gp$lik$likType, ") doesn't support Laplace-Fisher approximation. Try to set fisher.options = (sampling = TRUE).")
	}
}



# maybe delete this? :
d23_Fisher_dispatcher <- function (FUN, gp, f, data = gp$data, hyperpar, fisher = NULL)
{
	stopifnot(!is.null(fisher))
	stopifnot(FUN %in% c("d2", "d2_dhyp", "d3"))
	if (!is(gp[["lik"]], 'likTempPhased')) # !!!! ??? maybe not for the sampling???
		stop("For the Laplace-Fisher method, you need likTempPhased likelihood tempate, see ?gp.")
	if (!is.list(fisher))
		stop("parameter fisher has to be a list")
	if (!is.null(fisher[["sampling"]]) && fisher$sampling == TRUE) {
		return(do.call(paste0(FUN, "_Fisher_sampling"), list(gp, f, data, hyperpar, fisher)))
	} else {
		if (!is.character(gp$lik$likType))
			stop("For the Laplace-Fisher method, lik part of the likTempPhased template must be a character")
		if (gp$lik$likType == "bern")
			return(do.call(paste0(FUN, "_Fisher_bern"), list(gp, f, data, hyperpar, fisher)))
		else
			stop("gp(lik$likType = '", gp$lik$likType, ") doesn't support Laplace-Fisher approximation. Try to set fisher.options = (sampling = TRUE).")
	}
}

