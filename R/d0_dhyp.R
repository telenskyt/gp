
# derivative of d0 w.r.t. hyperparameter i
# i.e., the derivative of (the total joint) log likelihood w.r.t hyperparameter i (returns just one number)
# i - index within the optimized (non-fixed) hyperparameters
d0_dhyp <- function (gp, f, data = gp$data, hyperpar, i)
{
	i <- gpHyperparIdx(gp, i) - min(which(gp$hyperpar$component == ".lik")) + 1 # get the index of the likelihood hyperparameter within all likelihood hyperparameters (including the fixed ones)
	#numLikHyperpar <- sum(gp$hyperpar$component == ".lik") # number of all likelihood hyperparameters (including the fixed ones); we don't exclude the fixed ones for the RTMB tape
	
	if (gp$negLogLik.reindex2main && gp$GP_factor != "1") { 
		stopifnot(gpDataHasMainTable(data))
		# reindex from the GP_factor to main table in the closure wrapper, so that also the automatic differentiation
		# gives vector of the dimension of GP_factor!
		cmb <- function(func, data) function(par) {
			fact_idx <- paste0(gp$GP_factor, "_idx")
			par$f <- par$f[data[[1]][[fact_idx]]]
			func(data, par)
		}
	}
	else 
		cmb <- function(func, data) function(par) func(data, par)
			
	# cmb - the closure trick somewhere from the RTMB docs, to tie this to particular data and prevent unnecessary headaches!
	#mstart(id = "MakeTape", mem_precise = TRUE)
	par <- c(hyperpar[[".lik"]], list(f = numeric(gp$GP_size)))
	F <- MakeTape(cmb(gp$negLogLik, data), par)
	#mstop(id = "MakeTape")

	res <- drop(F$jacobian(unlist(par)))[i] # unlist() needed for the $jacobian() method
	if (!is.numeric(res) || !all(is.finite(res)))
		stop("d0_dhyp() of the user defined likelihood (negLogLik parameter to gp()) is not numeric: ", res)
	if (length(res) != 1)
		stop("length of the d0_dhyp() of user defined likelihood (negLogLik parameter to gp()) is ", length(res), ", should be 1")
	-res	
}
