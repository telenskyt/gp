#' Predict from a fitted Gaussian process model
#'
#' Predict from a fitted Gaussian process model.
#'
#' @param newdata object of class gpData, ***unscaled*** (i.e. not passed through \code{gpDataPrepare()}) data! They must be on the same scale as
#' the unscaled training dataset used for fitting the model. If \code{NULL} (the default), predictions will be made on the training dataset (which
#' is going to be super fast if conf.int = FALSE)
#' @param gp object of class gp, the fitted Gaussian process model
#' @param hyperpar (optional) lists of hyperparameter values, as returned by \code{gpHyperparList()}; default is to use the hyperparameters from the \code{gp} object.
#' @template param-components
#'
#' @param comp_missing how to compensate for missing components (if there are any), in mean and variance of the predictions:
##' \describe{
##' \item{\code{"none"}}{no compensation. In this case, the interpretation of the result's variance isn't much clear (!!!??? why?).}
##' \item{\code{"avg"}}{the missing components are replaced by "averaging over them" - more
##'		precisely, it predicts for the average effects of these components - for these, the variance/CI corresponds to confidence interval
##'		(uncertainty around the mean predictions) and not the prediction interval - see http://www.sthda.com/english/articles/40-regression-analysis/166-predict-in-r-model-predictions-and-confidence-intervals/
##'		For models which have grouping factors, be aware of how the average is weighted - see CAUTION CAU001 in the code.
##' }
#' }
#' Thus, the most reasonable setting now is \code{"avg"}; unless one cares only about the prediction mean (`f`) of the latent variable, in which case
#' a simple "none" would also do the job.
#'
#' @param w providing this optional argument will calculate average prediction(s) of all data points in this dataset, given these weights. Can
#' be either a vector of the same length as the size of \code{newdata} (in which case a single average is calculated), or,
#' if more different averages with different weights are desired, `w` can be a matrix where each column corresponds to one average - one set of weights
#' (i.e. number of rows of `w` must be equal to the size of \code{newdata}, number of columns equal to the number of different averages to be calculated).
#'
#' @param groupMeans optional, factor of the size of the predicted dataset. If provided, mean predictions for each factor level are calculated.
#'
#' @param type character; type of prediction:
#' \describe{
#' \item{\code{"latent"}}{predictions directly from the Gaussian Process itself (\code{f} and optionally \code{f_SE}), always in the dimension of the GP (\code{gp$GP_factor})}
#' \item{\code{"terms"}}{predictions passed through first stage of the likelihood template, that is, including GP (\code{f} and optionally \code{f_SE}) along with all the other terms used in the likelihood template. If \code{gp$lik.reindex2main} is \code{TRUE},
#'						the predictions will be reindexed to the dimension of the main table.}
#' \item{\code{"response"}}{predictions passed through the first two stages of the likelihood template, that is, after applying the inverse link function. If \code{gp$lik.reindex2main} is \code{TRUE},
#'						the predictions will be reindexed to the dimension of the main table.}
#' }
#'
#' @param conf.int logical; should confidence intervals be calculated? If \code{FALSE}, no confidence intervals will be calculated,
#' 	which makes the prediction much faster (but just in case both \code{se.fit} and \code{cov.fit} are also \code{FALSE}). 
#'	If \code{TRUE}, \code{se.fit} will be automatically set to \code{TRUE}, since the standard errors are needed for the CI calculation.
#' @param conf.level numeric; confidence level for the confidence intervals. Default is 0.95.
#' @param cov.fit logical; if \code{TRUE}, the covariance matrix of the latent \code{f} (the Gaussian process) will also be returned, in the original dimension of the GP. Note that the covariance matrix can be very memory consuming for large datasets.
#' Default is \code{FALSE}.
#' @param link character, optional. In case of \code{type = "response"}, what should be the link function to take inverse of for calculating the derived quantity (response scale)? Character as passed to \code{make.link()}.
#' @param parname character, optional. In case of \code{type = "response"}, what should be the name of the parameter which is the derived quantity on the response scale?
#' @param maxn maximum dataset size (along the \code{gp$GP_factor}) to fit at once; if the dataset is larger, it will be split into chunks of size \code{maxn}
#' before predictions are calculated. Use \code{Inf} to disable the splitting entirely. Choosing a suitable value may speed up the computation when \code{se.fit = TRUE} and \code{cov.fit = FALSE},
#' because the Kxx matrices are then smaller (that's probably the only part where splitting helps) - but testing is needed. (Note: splitting disabled for now
#' by default, since in some models we would need to specify along which factor to split - cannot always use gp$GP_factor)
#'
#' @param pred.sims (currently unused) number of simulations to be used for calculating the response scale predictions and their CIs. Default is 100000.
#' Higher values give more accurate results, but are slower. Ignored if \code{conf.int = FALSE}.
#' @param Kx.cache optional, object returned by \code{K_cache()} function, to speed up repeated calls to \code{predict.gp()} with the same \code{newdata} by
#' caching parts of the K(training_data, newdata) matrix.
#' @param Kxx.cache optional, object returned by \code{K_cache()} function, to speed up repeated calls to \code{predict.gp()} with the same \code{newdata} by
#' caching parts of the K(newdata, newdata) matrix.
#'
#' @returns A matrix of predictions, with columns \code{f} and \code{f_SE} (if \code{se.fit = TRUE}). If \code{cov.fit = TRUE}, returns a named list, \code{pred} will be
#'	the mentioned matrix and \code{cov} will be the full covariance matrix. For the dimension of the prediction, see the \code{type} parameter.
#'
#' @export

# pred.sims (timingy u plotu: = 1e6 => 24s; 1e5 => 12.58 ; 1e4 => 11s -> volim 1e5 jako default)
#	=> now is obsolete, we are using numerical integration. Whether or not it's used depends on conf.int.

# Note: weighting is done before predicting. Which is same as after, as it is all linear, but probably faster.

#
# groupMeans -

# ... - passed to model_expand_predictions()

predict.gp <- function(gp, newdata = NULL, type = c('latent', 'terms', 'response'),	se.fit = FALSE, cov.fit = FALSE, conf.int = FALSE, conf.level = 0.95,
						components = NULL, comp_missing = c("avg", "none"), w = NULL, groupMeans = NULL,
						hyperpar = gpHyperparList(gp),
						maxn = Inf, pred.sims = 100000,
						Kx.cache = NULL, Kxx.cache = NULL, ...)
{
	if (is.null(gp[["fit"]]))
		stop("Model object has not been fit yet: you need to call gpFit() first")
	comp_missing <- match.arg(comp_missing)

	need(gp, "data")
	need(gp, "obsdata") #!!! na co?

	components <- validate_components(gp, components) # model components (components of covariance matrix)
	all_components_used <- all(names(gp$covComp) %in% components)

	type <- match.arg(type)
	if (type != "latent" && (!is.null(w) || !is.null(groupMeans))) {
		warning("type != 'latent' together with w or groupMeans needs review; likelihood-template predictions may not be correct in this case.")
	}

	# do a backward compatibility check
	dots <- match.call(expand.dots = FALSE)$...
	if ("CI" %in% names(dots)) {
		stop("CI argument has been replaced by conf.int and conf.level")
	}

	# set up data
	same <- FALSE
	if (is.null(newdata)) {
		# use already set up inference data if not specified
		x_new <- gp$data
		if (all_components_used)
			same <- TRUE
	} else {
		# data must not be scaled!!!
		# use some hints to detect if user by mistake supplied scaled data
		if (gpDataIsScaled(newdata))
			stop("newdata must not be scaled for predict()! predict() will scale them. The newdata must be on the same scale as the unscaled training data.")
		# convert any ints to numerics (!!!2025: is it needed?)

		# at least check all if tables have the same columns as in the training dataset
        for (tbl_name in names(newdata)) {
            if (!tbl_name %in% names(gp$obsdata)) {
				stop("newdata contains table '", tbl_name, "' which was not present in the training data")
			} else {
				#if (!all(names(newdata[[tbl_name]]) %in% names(gp$obsdata[[tbl_name]]))) {
				#	stop("newdata$", tbl_name, " contains columns not present in training data")
				#}
				if (!all(names(gp$obsdata[[tbl_name]] %in% names(newdata[[tbl_name]])))) {
					stop("newdata$", tbl_name, " is missing columns present in training data")
				}
            }
        }

		x_new <- gpDataPrepare(gp, newdata, scale.as = gp$data) # it will scale the newdata in exactly the way that training data set was
	}

	# check confidence intervals
	if (!is.logical(conf.int) || length(conf.int) != 1L || is.na(conf.int)) {
		stop("conf.int must be TRUE or FALSE")
	}
	if (!is.numeric(conf.level) || length(conf.level) != 1L || is.na(conf.level) || conf.level >= 1 || conf.level <= 0) {
		stop("conf.level must be a number between 0 and 1")
	}
	if (conf.int) {
		se.fit <- TRUE
	}

	#if (is.null(maxn)) maxn <- ceiling(gpDataSize(x_new, gp$GP_factor)  / 10) # we are splitting the newdata, not the training dataset; ceiling needed instead of round here

	# first, put up together the latent prediction
	mstart(id = "whole predict", mem_precise = TRUE)
	mstart(id = "whole pred", mem_precise = TRUE)
	pred <- pred(gp, x_new, same = same, hyperpar = hyperpar, components = components, comp_missing = comp_missing,
		w = w, groupMeans = groupMeans, se.fit = se.fit, cov.fit = cov.fit, maxn = maxn, Kx.cache = Kx.cache, Kxx.cache = Kxx.cache
	)
	cat("pred() took ")
	mstop(id = "whole pred")
	pred.cov <- NULL
	if (cov.fit) {
		pred.cov <- pred$cov
		pred <- pred$pred
	}
	# pred is now always a matrix, with cols "f" and optionally also "f_SE"

	pred <- gpPredictFromLatent(gp = gp, pred = pred, data = x_new, type = type,
		hyperpar = hyperpar, se.fit = se.fit, cov.fit = cov.fit, pred.cov = pred.cov,
		conf.int = conf.int, conf.level = conf.level)

	cat("returning memory - gc() took ")
	mstart(id = "gc")
	gc() # Return the memory. Important! :-)
	mstop(id = "gc")
	cat("whole predict() took ")
	mstop(id = "whole predict")
	if (!cov.fit)
		stopifnot(is.matrix(pred))
	else
		stopifnot(is.list(pred), is.matrix(pred$pred))
	pred
}
