
#' Construct gaussian process model
#'
#' @param f formula for the covariance matrix
#' @param data object of class gpData, created by the gpData() function
#' @param logLik likelihood function taking two parameters: `data` (object of class gpData), and \code{f}, a numeric vector - result of the gaussian process.
#'			The \code{f} vector can be either of the dimension of the gaussian process, or the main table - see the \code{logLik.rescale} parameter.
#' @param logLik.rescale should the \code{f} parameter of the \code{logLik} function be rescaled to the dimension of the main table?
#'			If \code{FALSE}, the \code{f} parameter will be kept at the dimension of the Gaussian Process (as reported by \code{gpDetermineSize()}).
#'
#' @returns object of class gp - the gaussian process model. Some important slots:
#'\describe{
#'  \item{$hyperpar}{hyperparameter table (see XXX). User is invited to tweak these settings.}
#'  \item{$obsdata}{observed training data, in the original scale, as user supplied them}
#'  \item{$data}{training data prepared for model fit - converted to matrices, and scaled (standardized) where appropriate}
#'}
#' @export
gp <- function(f, data, logLik, logLik.rescale = TRUE)
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
	gp$ll.rescale_to_main <- FALSE
	if (logLik.rescale && gp$GP_factor != "1") { # rescaling is desired to take place
		if (!gpDataHasMainTable(gp$data))
			stop("The parameter logLik.rescale = TRUE, but I cannot rescale from the factor ", gp$GP_factor, ", at which gaussian process is running, to the main table, since it is missing in the training data. Supply the main table or consider setting logLik.rescale = FALSE.")
		gp$ll.rescale_to_main <- TRUE
	}

	gp
}
