
# bernoulli likelihood; invoked when likelihood is given as string - "bern"
likBern <- function(par, log = FALSE, sum = FALSE)
{	
	getAll(par, warn=FALSE)
	stopifnot(!sum || log) # sum => log
	if (is.matrix(y) && !is.matrix(p) && isTRUE(log) && isTRUE(sum)) {
		# dbinom() would work perfectly anyway for x being matrix and prob vector,
		# but this little optimization - pooling the likelihood across rows - significantly speeds up any RTMB::MakeTape created from this vectorized case! (e.g. in Fisher_sampling)		
		# but it requires that the function calculates the sum already (sum = TRUE)
		y.rs <- rowSums(y)
		stopifnot(length(y.rs) == length(p))
		
		sum(y.rs*log(p) + (ncol(y)-y.rs)*log1p(-p))
			# the only problem here is for p exactly 0 or 1, the result will be 0*log(0) = NaN, see https://chatgpt.com/c/6a19f4be-1e6c-83eb-a785-3c2c2f0bae1d
			# dbinom() has this handled
			# but lets ignore it for now, hopefully it won't happen
	} else {
		res <- dbinom(x = y, size = 1, prob = p, log = log)
		if (isTRUE(log) && isTRUE(sum))
			sum(res)
		else
			res
	}
}

# samples ~ columns
# had to write my own; using RTMB's $simulate() was sketchy (the y <- OBS(y) even converted matrix to vector) and slow
simBern <- function (par, samples)
{
	list(y = matrix(rbinom(n = samples*length(par$p), size = 1, prob = par$p), ncol = samples))
}



#' Evaluation of predictive performance for bernoulli (0/1) likelihood
#'@param par the output of processing stage (stage 3) of the likelihood template. For the bernoulli likelihood in particular, it is 
#' a list with two elements:
#'\describe{
#'\item{p}{numeric vector - predicted probability}
#'\item{y}{bernoulli response variable (vector of 0s and 1s)}
#'}
#'@param par.null the same as \code{par}, but for the null model (with intercept only). Only the \code{p} will be used from this list.
#'
#'@returns list with the calculated metrics:
#'	- AUC - scale from 0 to 1, random prediction is 0.5
#'
#'	- TSS - True skill statistics - scale from 0 to 1, random prediction is 0 \insertCite{allouche_assessing_2006}{gp}
#'
#'	- R2_dev - pseudo R-squared based on binomial deviance. Scale from 0 to 1, 0 is null model (with intercept only)
#'
#'  - R2_LR - McFadden's pseudo R-squared based on likelihood ratio. Scale from 0 to 1, 0 is null model (with intercept only) \insertCite{office_of_advanced_research_computing_university_of_california_los_angeles_faq_2011}{gp}
#'
#'	- R2_LRNCU - Nagelkerke, Cragg & Uhler's pseudo R-squared based on likelihood ratio, corrected for sample size. Scale from 0 to 1, 0 is null model (with intercept only) \insertCite{office_of_advanced_research_computing_university_of_california_los_angeles_faq_2011}{gp}
#'
#' Perfect prediction (the limit) has always value 1 in all these statistics.
#'
#' Note that the AUC and TSS scale relative to random prediction, while the pseudo R-squared statistics scale relative to a null model (with intercept only).
#' Since the null model has already more information than the completely random prediction, it can scale with relatively high AUC and TSS, while showing up as 0 in the R2 statistics.
#' This is the reason why the AUC and TSS statistics will always "look better" than the R2_* statistics.
#'
#' Note that for models with more complicated likelihood the AUC and TSS might not make sense. The R2_* statistics are more general, and thus more robust sense-wise.
#'
#'@importFrom WeightedROC WeightedAUC WeightedROC
#'@importFrom ecospat ecospat.max.tss
#'@importFrom Rdpack reprompt
#'@export
#'@references
#'	\insertAllCited{}
#
# Notes: 2026-01: 
#	- I see that the p.null is necessary, in the case of Lapwing nest survival, even in the null model, the `p` probability will need to have different number of days based on length
# 	  of the time period the nest failed!
#	- I don't provide Cox & Snell variant, since it is not correctly scaled (the NCU variant is the correction)
predPerfBern <- function(par, par.null)
{
	p <- par$p
	y <- par$y
	p.null <- par.null$p
	
	AUC <- WeightedAUC(WeightedROC(p, y)) # same as in the Czech Atlas (Stastny et al., 2021)
	TSS <- ecospat.max.tss(p, y)$max.TSS

	if (is.null(p.null)) {
		p.null <- sum(y)/length(y) # for TSS, it is searching for optimal threshold; here, to get the "null deviance", we just estimate it!!
	}
	stopifnot(length(y) == length(p))
	stopifnot(length(y) == length(p.null) || length(p.null) == 1)
	
	dev <- devianceBin(y, 1, p)
	nullDev <- devianceBin(y, 1, p.null)
	
	NLL <- -sum(dbinom(x = y, size = 1, prob = p, log = TRUE))
	nullNLL <- -sum(dbinom(x = y, size = 1, prob = p.null, log = TRUE))
	N <- length(y)
	
	
	list(
		AUC = AUC, 
		TSS = TSS, 
		R2_dev = 1 - dev / nullDev, # same as the R2_dev in the Czech Atlas (Stastny et al., 2021)
		R2_LR = 1 - NLL/nullNLL, # McFadden's, see also https://stats.idre.ucla.edu/other/mult-pkg/faq/general/faq-what-are-pseudo-r-squareds/ and also package rsq
		R2_LRNCU = (1 - exp(-2/N*(nullNLL - NLL)))/(1 - exp(-2/N*nullNLL)), # Nagelkerke / Cragg & Uhler’s, see also https://stats.idre.ucla.edu/other/mult-pkg/faq/general/faq-what-are-pseudo-r-squareds/ and also package rsq.
			# likelihood ratios:
			#	0 = null model (with intercept only)
			#	1 = perfect prediction 
		NLL = NLL,
		nullNLL = nullNLL
	)
}

# deviance of binomial glm!
# https://www.physicsforums.com/threads/deviance-of-binomial-generalized-linear-model.430009/
devianceBin <- function (y, n, p) {
	small_number <- 1e-13 # nastesti moc nezalezi na tom, jak maly to cislo bude :)
	yhat <- n*p
	#2*sum(y*log((y+small_number)/yhat)+(n-y)*log((n-y+small_number)/(n-yhat)))
	2*sum(y*log((y+small_number)/(yhat+small_number))+(n-y)*log((n-y+small_number)/(n-yhat+small_number)))
}