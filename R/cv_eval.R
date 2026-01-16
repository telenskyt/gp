

#' Evaluation of cross-validation statistics for bernoulli (0/1) response variable and predicted probability
#' 
#' Calculates AUC, TSS, and R2_dev (binomial deviance)
#'@param p numeric vector - predicted probability
#'@param y bernoulli response variable (vector of 0s and 1s)
#'@returns list with the calculated metrics
#'
#'@importFrom WeightedROC WeightedAUC WeightedROC
#'@importFrom ecospat ecospat.max.tss
#'@export
cv_eval_bern <- function(p, y)
{
	AUC <- WeightedAUC(WeightedROC(p, y)) # same as in the Czech Atlas (Stastny et al., 2021)
	TSS <- ecospat.max.tss(p, y)$max.TSS

	p_Null <- sum(y)/length(y) # for TSS, it is searching for optimal threshold; here, to get the "null deviance", we just estimate it!!
	dev <- devianceBin(y, 1, p)
	nullDev <- devianceBin(y, 1, p_Null)
	R2_dev <- 1 - dev / nullDev # same as the R2_dev in the Czech Atlas (Stastny et al., 2021)
	
	list(AUC = AUC, TSS = TSS, R2_dev = R2_dev)
}


# deviance of binomial glm!
# https://www.physicsforums.com/threads/deviance-of-binomial-generalized-linear-model.430009/
devianceBin <- function (y, n, p) {
	small_number <- 1e-13 # nastesti moc nezalezi na tom, jak maly to cislo bude :)
	yhat <- n*p
	#2*sum(y*log((y+small_number)/yhat)+(n-y)*log((n-y+small_number)/(n-yhat)))
	2*sum(y*log((y+small_number)/(yhat+small_number))+(n-y)*log((n-y+small_number)/(n-yhat+small_number)))
}