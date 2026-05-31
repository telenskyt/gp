

# Fisher template for any distribution, using sampling!
#	- modification of neg. log lik. template
#	- hessian of this function will give the Fisher information matrix :)
#
# This version is now super fast!

Fisher_sampling <- function (gp, f, data, hyperpar, fisher)
{
	samples <- 1000
	if (!is.null(fisher[["samples"]]))
		samples <- fisher$samples
	# tento sampling funguje! Ale bude mozna pomalej
	par <- gpGetParForLikTemplate(gp, f, data, hyperpar)
	par2 <- gp$lik$stages(data, par, stages = 1:3)

	simdat <- gp$lik$simulate(par2, samples)
		
	fisher_fun <- function (data, par) {
		par2 <- gp$lik$stages(data, par, stages = 1:3)
			# !!! tady by slo urcite jeste zrychlit vektorizovanim teto faze
		
		par3 <- modifyList(par2, simdat)
		
		nllsum <- gp$lik$stages(data, par3, stages = 4:5)
		return(nllsum/samples)
	}
	attr(fisher_fun, "eliminate") <- TRUE
	attr(fisher_fun, "optimize") <- TRUE
	fisher_fun
}
