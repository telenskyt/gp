
#' Construct gaussian process model
#'
#' @param f formula for the covariance matrix
#' @param data object of class gpData, created by the gpData() function
#' @param negLogLik user-defined negative log likelihood function, taking two parameters: `data` (object of class gpData), and \code{par}, a list of parameters (numeric scalars or vectors).
#'		The \code{par} parameters list always includes \code{f}, a numeric vector - result of the gaussian process. The \code{par$f} vector can be either of the dimension of the gaussian process, 
#'		or the main table - see the \code{negLogLik.reindex2main} parameter. The \code{par} parameters list can optionally include other parameters, see the \code{negLogLik.hyperpar} argument below.
#' @param negLogLik.hyperpar optional, a named list of extra hyperparameters for the \code{negLogLik} function (named list of numeric scalars or vectors); this doesn't include the special 
#' 		parameter \code{f} (the Gaussian process).
#' @param negLogLik.formula optional, formula related to the supplied \code{negLogLik} function. This is not used in any way, except it is stored in the object and printed in the summary,
#'		so it has only informational purpose.
#' @param negLogLik.reindex2main should the \code{par$f} parameter of the \code{negLogLik} function be reindexed to the rows of the main table?
#'			If \code{FALSE}, the \code{par$f} parameter will be kept at the dimension of the Gaussian Process (as reported by \code{gpDetermineSize()}).
#'			If the Gaussian process is running at the dimension of the main table, then it does not matter.
#'
#' @returns object of class gp - the gaussian process model. Some important slots:
#'\describe{
#'  \item{$hyperpar}{hyperparameter table (see XXX). User is invited to tweak these settings.}
#'  \item{$obsdata}{observed training data, in the original scale, as user supplied them}
#'  \item{$data}{training data prepared for model fit - converted to matrices, and scaled (standardized) where appropriate}
#'}
#' @export
gp <- function(f, data, negLogLik, negLogLik.hyperpar = NULL, negLogLik.formula = NULL, negLogLik.reindex2main = TRUE,
	predictor.fun = NULL, response.parname = NULL, link = NULL)
{
	gp <- gpFormula(f)
	gp$covFormula <- f
	gp$predictor.fun <- predictor.fun
	gp$response.parname <- response.parname
	gp$link <- link

	stopifnot(class(data) == "gpData")

	gp$obsdata <- data

	stopifnot(gpDataCheckReq(gp, data))

	gp$data <- gpDataPrepare(gp, data)

	gp$covComp_df <- gpComponentsTable(gp)

	gp_size <- gpDetermineSize(gp)
	gp$GP_size <- gp_size$size
	gp$GP_factor <- gp_size$fact # bude vzdy character string ruzny od "", NA, NULL, viz gpSize()

	stopifnot(is.function(negLogLik))
	gp$negLogLik <- negLogLik
	if (is.null(negLogLik.hyperpar))
		negLogLik.hyperpar <- list()
	gp$negLogLik.hyperpar <- negLogLik.hyperpar
	if ("f" %in% names(gp$negLogLik.hyperpar))
		stop("The negLogLik.hyperpar cannot contain parameter named 'f' - this is reserved for the Gaussian process.")
	gp$negLogLik.formula <- negLogLik.formula
	gp$negLogLik.reindex2main <- negLogLik.reindex2main
	if (gp$negLogLik.reindex2main && gp$GP_factor != "1") { # the reindexing is desired to take place
		if (!gpDataHasMainTable(gp$data))
			stop("The parameter negLogLik.reindex2main = TRUE, but I cannot reindex from the factor ", gp$GP_factor, ", at which gaussian process is running, to the main table, since it is missing in the training data. Supply the main table or consider setting negLogLik.reindex2main = FALSE.")
	}

	gp$hyperpar <- gpHyperparDefaults(gp)
	gpHyperparCheckAll(gp)

	gp
}
