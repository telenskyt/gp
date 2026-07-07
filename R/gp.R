
#' Construct gaussian process model
#'
#' @param f formula for the covariance matrix
#' @param data object of class gpData, created by the gpData() function
#' @param lik an object of class likTemp, representing the user-defined template of the likelihood function \code{p(y|f)}.
#' 
#' ORIG DOC TO BE DELETED - MAYBE USE PART OF THESE FOR THE DOCUMENTATION OF CLASS likTemp
#' either a single function, or a named list of three functions (three-stage template); these are two ways to supply 
#' the user-defined negative log likelihood function p(y|f):
#' \enumerate{
#'	\item \code{negLogLik} is a single function \code{f(data, par)}, taking two parameters: 
#'	\describe{
#'		\item{data}{object of class \code{gpData}, containing both predictors and a response variable}
#'		\item{\code{par}}{a list of parameters (numeric scalars or vectors).
#'			Always includes \code{f}, a numeric vector - result of the gaussian process. The \code{par$f} vector can be either of the dimension of the gaussian process, 
#'			or the main table - see the \code{negLogLik.reindex2main} parameter. The \code{par} parameters list can optionally include other parameters, see the \code{negLogLik.hyperpar} argument below.
#'   	}
#'	}
#'	\item Logic is similar as in the previous case, but the \code{negLogLik} is supplied as a chain of three consecutive functions (three-stage template). Then, \code{negLogLik} is a named list of three functions:
#'	\describe{
#'		\item{\code{pred.link}}{either character or \code{function(data, par)}. This is predictor & link component. Processes the parameters into other parameters; 
#'			is also used by \code{predict(type = "response")} for predictions,
#'			therefore, the \code{data} might not contain the response variable (thus, the response variable must not be used by this function). Should be either a character string representing 
#'			the name of a link function	(this will be the most frequent case; see !!! for supported link functions), or a \code{function(data, par)} (see above for parameter description !!!).
#'			This function will return a named list of parameters \code{par} that will be passed on to the next stage.}
#'		\item{\code{process}}{\code{function(data, par)}, this function processes the parameters and the \code{data} (which always contain the response variable here) and shapes both
#'			parameters and data for particular likelihood function. This function should return a named list required by particular likelihood function (see !!!).
#'			This component is used by cross-validation and Laplace-Fisher approximation.}
#'		\item{\code{lik}}{function evaluating the likelihood. Either character representing name of likelihood function (see !!! for supported values), or custom \code{function(x)} taking \code{x} from 
#'		the previous stage.}
#' }
#' }
#' @param lik.hyperpar optional, a named list of extra likelihood hyperparameters for the \code{negLogLik} function (named list of numeric scalars or vectors); this doesn't include the special 
#' 		parameter \code{f} (the Gaussian process). This effectively defines the likelihood hyperparameters, named ".lik" in the hyperparameter table.
#' @param lik.formula optional, formula related to the supplied \code{negLogLik} function. This is not used in any way, except it is stored in the object and printed in the summary,
#'		so it has only informational purpose.
#' @param lik.reindex2main should the \code{par$f} parameter of the \code{negLogLik} function be reindexed to the rows of the main table?
#'			If \code{FALSE}, the \code{par$f} parameter will be kept at the dimension of the Gaussian Process (as reported by \code{gpDetermineSize()}).
#'			If the Gaussian process is running at the dimension of the main table, then it does not matter.
#'
#' @param W.type the type of the hessian matrix (\code{W}) of the negative log likelihood function (as defined by the \code{lik} parameter)
#' \describe{
#' \item{\code{"diag"}}{(default) - the hessian matrix (\code{W}) is diagonal. Most common case, fastest computation. This is usually the case when the user-defined likelihood factorizes over individual \eqn{f_i}'s.}
#' \item{\code{"bdiag"}}{the hessian matrix (\code{W}) is block-diagonal. This is usually the case when the user-defined likelihood factorizes over mutually disjoint sets of \eqn{f_i}'s.}
#' \item{\code{"other"}}{the hessian matrix (\code{W}) is a general matrix.}
#' }
#' One should keep the setting as specific as possible to keep the computation efficient.
#'
#' @returns object of class gp - the gaussian process model. Some important slots:
#'\describe{
#'  \item{$hyperpar}{hyperparameter table (see XXX). User is invited to tweak these settings.}
#'  \item{$obsdata}{observed training data, in the original scale, as user supplied them}
#'  \item{$data}{training data prepared for model fit - converted to matrices, and scaled (standardized) where appropriate}
#'}
#' @export
gp <- function(f, data, lik, lik.hyperpar = NULL, lik.formula = NULL, lik.reindex2main = TRUE,
	W.type = c("diag", "bdiag", "other"),
	predictor.fun = NULL, response.parname = NULL, link = NULL)
{
	gp <- gpFormula(f)
	gp$covFormula <- f
	gp$predictor.fun <- predictor.fun
	gp$response.parname <- response.parname
	gp$link <- link

	gp <- gpSetData(gp, data)

	stopifnot(is(lik, "likTempPhased"))
	gp$lik <- lik
	gp$lik.class <- lik$getRefClass()$className
	gp$lik.constr.args <- lik$constr.args
	if (is.null(lik.hyperpar))
		lik.hyperpar <- list()
	gp$lik.hyperpar <- lik.hyperpar
	if ("f" %in% names(gp$lik.hyperpar))
		stop("The lik.hyperpar cannot contain parameter named 'f' - this is reserved for the Gaussian process.")
	gp$lik.formula <- lik.formula
	gp$lik.reindex2main <- lik.reindex2main
	if (gp$lik.reindex2main && gp$GP_factor != "1") { # the reindexing is desired to take place
		if (!gpDataHasMainTable(gp$data))
			stop("The parameter lik.reindex2main = TRUE, but I cannot reindex from the factor ", gp$GP_factor, ", at which gaussian process is running, to the main table, since it is missing in the training data. Supply the main table or consider setting lik.reindex2main = FALSE.")
	}

	gp$W.type <- match.arg(W.type)
	gp$hyperpar <- gpHyperparDefaults(gp)
	gpHyperparCheckAll(gp)
	gp$.pkg.version <- pkg_build_info()

	gp
}

# function to set a data as a training dataset to a gp object
#
# used in gp(), but also in gpFitCV() - is this a sign that the interface 
# should have been different, not passing data to gp() but to gpFit() and gpFitCV()?

gpSetData <- function (gp, data)
{
	stopifnot(class(data) == "gpData")

	gp$obsdata <- data

	stopifnot(gpDataCheckReq(gp, data))

	gp$data <- gpDataPrepare(gp, data)

	gp$covComp_df <- gpComponentsTable(gp)

	gp_size <- gpDetermineSize(gp)
	gp$GP_size <- gp_size$size
	if (!is.null(gp[["GP_factor"]]))
		stopifnot(gp$GP_factor == gp_size$fact) # should still be the same
		
	gp$GP_factor <- gp_size$fact # bude vzdy character string ruzny od "", NA, NULL, viz gpDetermineSize()

	gp
}

