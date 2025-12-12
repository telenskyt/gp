
#' Construct gaussian process model
#'
#' @param f formula for the covariance matrix
#' @param data object of class gpData, created by the gpData() function
#' @param logLik likelihood function taking two parameters: `data` (object of class gpData), and \code{f}, a numeric vector - result of the gaussian process.
#'			The \code{f} vector can be either of the dimension of the gaussian process, or the main table - see the \code{logLik.reindex2main} parameter.
#' @param logLik.reindex2main should the \code{f} parameter of the \code{logLik} function be reindexed to the rows of the main table?
#'			If \code{FALSE}, the \code{f} parameter will be kept at the dimension of the Gaussian Process (as reported by \code{gpDetermineSize()}).
#'			If the Gaussian process is running at the dimension of the main table, then it does not matter.
#'
#' @returns object of class gp - the gaussian process model. Some important slots:
#'\describe{
#'  \item{$hyperpar}{hyperparameter table (see XXX). User is invited to tweak these settings.}
#'  \item{$obsdata}{observed training data, in the original scale, as user supplied them}
#'  \item{$data}{training data prepared for model fit - converted to matrices, and scaled (standardized) where appropriate}
#'}
#' @export
gp <- function(f, data, logLik, logLik.reindex2main = TRUE)
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
	gp$logLik.reindex2main <- llogLik.reindex2main
	if (gp$logLik.reindex2main && gp$GP_factor != "1") { # the reindexing is desired to take place
		if (!gpDataHasMainTable(gp$data))
			stop("The parameter logLik.reindex2main = TRUE, but I cannot reindex from the factor ", gp$GP_factor, ", at which gaussian process is running, to the main table, since it is missing in the training data. Supply the main table or consider setting logLik.reindex2main = FALSE.")
	}

	gp
}
