# derivative of d1 w.r.t. hyperparameter i
# i.e. second-order partial derivative of log likelihood given f w.r.t f w.r.t hyperpar i (returns vector)
# i - index within the optimized (non-fixed) hyperparameters, must be a likelihood hyperparameter
d1_dhyp <- function (gp, f, data = gp$data, hyperpar, i)
{
	# validation of i
	stopifnot(gp$hyperpar$component[gpHyperparIdx(gp, i)] == ".lik")
	stopifnot(all(gp$hyperpar$component[min(which(gp$hyperpar$component == ".lik")):gpHyperparIdx(gp, i)] == ".lik"))
	i <- gpHyperparIdx(gp, i) - min(which(gp$hyperpar$component == ".lik")) + 1 # get the index of the likelihood hyperparameter within all likelihood hyperparameters (including the fixed ones)
	numLikHyperpar <- sum(gp$hyperpar$component == ".lik") # number of all likelihood hyperparameters (including the fixed ones); we don't exclude the fixed ones for the RTMB tape
	stopifnot(i >= 1)
	stopifnot(i <= numLikHyperpar)
	
	if (gp$lik.reindex2main && gp$GP_factor != "1") { 
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
	F <- MakeTape(cmb(gp$lik$nll, data), par)
	#mstop(id = "MakeTape")
	
	F1 <- F$jacfun()

	F1T <- MakeTape(function (f) { 
		"[<-" <- ADoverload("[<-") 
		"c" <- ADoverload("c")		
		#par[["f"]] <- f
		#F$jacobian(unlist(par))[i] 
		
		#F1(c(hyperpar[[".lik"]], list(f = f)))[i] # toto hazelo chybu
		F1(c(unlist(hyperpar[[".lik"]]), f))[i] # tohle je reseni te chyby, pomohl github copilot, ale nechapu proc tady najednou je potreba unlist, to nikdy nebylo!
	}, f) # tape to take just the hessian diagonal (and just for the f)

	F1T$simplify("eliminate")
	F1T$simplify("optimize")
	#F1T$atomic() # doesn't speed it up

	F2 <- F1T$jacfun()	
	res <- drop(F2(f))
	
	if (!is.numeric(res) || !all(is.finite(res)))
		stop("second derivative of the user defined likelihood (negLogLik parameter to gp()) is not numeric: ", res)
	if (length(res) != gp$GP_size)
		stop("length of the second derivative the of user defined likelihood (negLogLik parameter to gp()) is ", length(res), ", should be ", gp$GP_size)
	-res
}
