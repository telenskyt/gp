
#' Construct gaussian process model
#'
#' @param f formula for the covariance matrix
#' @param data object of class gpData, created by the gpData() function
#' @param logLik blikelihood function taking two parameters: `data` (object of class gpData), and `f`, a numeric vector - result of the gaussian process
#'
#' @returns object of class gp - the gaussian process model. Some important slots:
#'\describe{
#'  \item{$hyperpar}{hyperparameter table (see XXX). User is invited to tweak these settings.}
#'  \item{$obsdata}{observed training data, in the original scale, as user supplied them}
#'  \item{$data}{training data prepared for model fit - converted to matrices, and scaled (standardized) where appropriate}
#'}
#' @export
gp <- function(f, data, logLik)
{
	gp <- gpFormula(f)
	gp$covFormula <- f

	stopifnot(class(data) == "gpData")

	gp$obsdata <- data

	stopifnot(gpDataCheckReq(gp, data))

	gp$data <- gpDataPrepare(gp, data)

	gp$covComp_df <- gpComponentsTable(gp)

	gp_size <- gpDetermineSize(gp)
	gp$GP_size <- gp_size$size
	gp$GP_factor <- gp_size$fact # bude vzdy character string ruzny od "", NA, NULL, viz gpSize()

	gp$hyperpar <- gpHyperparDefaults(gp)

	stopifnot(is.function(logLik))
	gp$ll <- logLik


	gp
}
