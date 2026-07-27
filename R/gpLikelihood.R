
# The functions d[0-3]* are now in separate files
#
# obecna note k temto fcim - zatim nechavam data jako separatni parametr, i kdyz by mozna stacilo brat proste data z gp$data
# (ale co cross-validace treba?)

gpGetParForLikTemplate <- function(gp, f, data, hyperpar)
{
	if (gp$lik.reindex2main && gp$GP_factor != "1") { 
		stopifnot(gpDataHasMainTable(data))	
		# reindex from the GP_factor to main table, as requested by the template
		# this should be called within the MakeTape template, so that also the automatic differentiation
		# gives vector of the dimension of GP_factor!
		fact_idx <- paste0(gp$GP_factor, "_idx")
		f <- f[data[[1]][[fact_idx]]]
	}
	par <- c(hyperpar[[".lik"]], list(f = f))
	par
}

#' Run the GP likelihood template
#'
#' Convenience function to call the GP likelihood template with the gp object defaults
#' (fitted latent values and the model's likelihood hyperparameters).
#' Probably will not be part of a normal workflow; servers the user who wants to "dig" deeper, maybe for some diagnostic purposes.
#'
#' @param gp GP model object
#' @param par Input for the first requested likelihood-template stage. If
#'	`NULL`, it is constructed from `gp$fit$f` and the `.lik` component of
#'	`hyperpar`, in the format expected by stage 1. Therefore, `par` may be
#'	`NULL` only when stage 1 is selected.
#' @param stages Integer vector selecting the stages to run. Available stages
#'	are 1 (terms), 2 (link), 3 (likelihood-specific processing),
#'	4 (element-wise likelihood), and 5 (negative sum of log likelihoods).
#'	The selected stages are run in this pipeline order; skipped stages are
#'	not evaluated.
#' @param hyperpar A list of hyperparameter values. Only its `.lik` component is used, and only when
#'	`par` is `NULL`.
#'
#' @details
#' When stage 1 is omitted, `par` must be supplied in the format expected by
#' the first selected stage, normally the output returned by the preceding
#' stage.
#'
#' @return The output of the last selected stage. Its type depends on that
#' stage and on the likelihood template. With the defaults, this is a numeric
#' scalar containing the negative log likelihood.
#'
#' @examples
#' \dontrun{
#' # gp_fit is a fitted gp model
#' gpLik(gp_fit)
#' gpLik(gp_fit, stages = 1:2)
#' }
#'
#' @seealso [gpHyperparList()], [gpPredictFromLatent()]
#' @export
gpLik <- function (gp, par = NULL, stages = 1:5, hyperpar = gpHyperparList(gp))
{
	data <- gp$data 
	if (is.null(par)) {
		if (!1 %in% stages)
			stop("if stage 1 is missing, par must be provided")
		if (is.null(gp[["fit"]]))
			stop("if gp$fit doesn't exist, `par` must be provided")
		stopifnot(!is.null(gp$fit[["f"]]))
		par <- gpGetParForLikTemplate(gp, f = gp$fit$f, data = data, hyperpar = hyperpar)
	}
		
	gp$lik$stages(data = data, par = par, stages = stages)
}
