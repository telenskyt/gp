

#' Evaluation of cross-validation statistics for bernoulli (0/1) response variable and predicted probability
#'
#'@param p numeric vector - predicted probability
#'@param y bernoulli response variable (vector of 0s and 1s)
#'@param p.null the probability predicted by a null model (with intercept only)
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
cv_eval_bern <- function(p, y, p.null)
{
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