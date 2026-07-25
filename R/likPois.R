
# poisson likelihood; invoked when likelihood is given as string - "pois"
likPois <- function(par, log = FALSE, sum = FALSE)
{	
	getAll(par, warn=FALSE)
	stopifnot(!sum || log) # sum => log
	if (is.matrix(y) && !is.matrix(lambda) && isTRUE(log) && isTRUE(sum)) {
		# dpois() would work perfectly anyway for x being matrix and lambda vector,
		# but this little optimization - pooling the likelihood across rows - significantly speeds up any RTMB::MakeTape created from this vectorized case!		
		# but it requires that the function calculates the sum already (sum = TRUE)
		y.rs <- rowSums(y)
		stopifnot(length(y.rs) == length(lambda))
		
		sum(y.rs*log(lambda) - ncol(y)*lambda) - sum(lgamma(y+1))
			# the only problem here is for lambda exactly 0, the result will be 0*log(0) = NaN
			# dpois() has this exact lambda == 0 boundary case handled, but it cannot recover eta if exp(eta) underflowed to 0
			# the numerically stable solution of this would be to pass the linear predictor directly (as `eta`) and use y.rs*eta - ncol(y)*lambda, 
			# but that would be complicated now given the current template design
			# so lets ignore it for now, hopefully it won't happen
	} else {
		res <- dpois(x = y, lambda = lambda, log = log)
		if (isTRUE(log) && isTRUE(sum))
			sum(res)
		else
			res
	}
}

# samples ~ columns
# had to write my own; using RTMB's $simulate() was sketchy (the y <- OBS(y) even converted matrix to vector) and slow
simPois <- function (par, samples)
{
	list(y = matrix(rpois(n = samples*length(par$lambda), lambda = par$lambda), ncol = samples))
}



#' Evaluation of predictive performance for poisson likelihood
#'@param par the output of processing stage (stage 3) of the likelihood template. For the poisson likelihood in particular, it is 
#' a list with two elements:
#'\describe{
#'\item{lambda}{numeric vector - predicted mean}
#'\item{y}{poisson response variable (counts)}
#'}
#'@param par.null the same as \code{par}, but for the null model (with intercept only). Only the \code{lambda} will be used from this list.
#'
#'@returns list with the calculated metrics:
#'	- R2 - squared ordinary correlation between observed and predicted values. Scale from 0 to 1
#'
#'	- R2_dev - pseudo R-squared based on poisson deviance. Scale from 0 to 1, 0 is null model (with intercept only)
#'
#'  - R2_LR - McFadden's pseudo R-squared based on likelihood ratio. Scale from 0 to 1, 0 is null model (with intercept only)
#'
#'	- R2_LRNCU - Nagelkerke, Cragg & Uhler's pseudo R-squared based on likelihood ratio, corrected for sample size. Scale from 0 to 1, 0 is null model (with intercept only)
#'
#' Perfect prediction (the limit) has always value 1 in all these statistics.
#'
#' Note that for models with more complicated likelihood the R2 might not make sense. The R2_* statistics are more general, and thus more robust sense-wise.
#'
#'@export
predPerfPois <- function(par, par.null)
{
	lambda <- par$lambda
	y <- par$y
	lambda.null <- par.null$lambda
	
	R2 <- cor(y, lambda)^2
	
	if (is.null(lambda.null)) {
		lambda.null <- mean(y) # here, to get the "null deviance", we just estimate it!!
	}
	stopifnot(length(y) == length(lambda))
	stopifnot(length(y) == length(lambda.null) || length(lambda.null) == 1)
	
	dev <- deviancePois(y, lambda)
	nullDev <- deviancePois(y, lambda.null)
	
	NLL <- -sum(dpois(x = y, lambda = lambda, log = TRUE))
	nullNLL <- -sum(dpois(x = y, lambda = lambda.null, log = TRUE))
	N <- length(y)
	
	
	list(
		R2 = R2, 
		R2_dev = 1 - dev / nullDev,
		R2_LR = 1 - NLL/nullNLL, # McFadden's, see also https://stats.idre.ucla.edu/other/mult-pkg/faq/general/faq-what-are-pseudo-r-squareds/ and also package rsq
		R2_LRNCU = (1 - exp(-2/N*(nullNLL - NLL)))/(1 - exp(-2/N*nullNLL)), # Nagelkerke / Cragg & Uhler's, see also https://stats.idre.ucla.edu/other/mult-pkg/faq/general/faq-what-are-pseudo-r-squareds/ and also package rsq.
			# likelihood ratios:
			#	0 = null model (with intercept only)
			#	1 = perfect prediction 
		NLL = NLL,
		nullNLL = nullNLL
	)
}

# deviance of poisson glm!
# using the exact same formula as in R - poisson()$dev.resids
deviancePois <- function (y, lambda) {
	lambda <- lambda + y*0
	res <- lambda
	pos <- y > 0
	res[pos] <- y[pos]*log(y[pos]/lambda[pos]) - (y[pos]-lambda[pos])
	2*sum(res)
}
