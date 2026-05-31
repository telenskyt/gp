
# Fisher template for bernoulli distribution
#	- modification of neg. log lik. template
#	- hessian of this function will give the Fisher information matrix :)

Fisher_bern <- function (gp, f, data = gp$data, hyperpar, fisher = NULL)
{
	par <- gpGetParForLikTemplate(gp, f, data, hyperpar)
	proc <- gp$lik$stages(data, par, stages = 1:3)
	const_p_equals_1 <- proc$p
	
	fisher_fun <- function(data, par) {
		proc <- gp$lik$stages(data, par, stages = 1:3) # evaluation of the stages terms, link & process
		-sum((1 - const_p_equals_1)*log1p(-proc$p) + const_p_equals_1 * log(proc$p))
			# the only problem here is for p exactly 0 or 1, the result will be 0*log(0) = NaN, see https://chatgpt.com/c/6a19f4be-1e6c-83eb-a785-3c2c2f0bae1d
			# dbinom() has this handled
			# but lets ignore it for now, hopefully it won't happen			
	}
	
	fisher_fun
}
