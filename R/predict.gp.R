#' Predict from a fitted Gaussian process model
#' 
#' Predict from a fitted Gaussian process model. 
#' 
#' @param newdata object of class gpData, ***unscaled*** (i.e. not passed through \code{gpDataPrepare()}) data! They must be on the same scale as 
#' the unscaled training dataset used for fitting the model. If \code{NULL} (the default), predictions will be made on the training dataset (which
#' is going to be super fast if CI = NULL)
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
#' \item{\code{"terms"}}{predictions including GP (\code{f} and optionally \code{f_SE}) along with all the other terms used in the likelihood template. If gp$lik.reindex2main = TRUE, 
#'						the predictions will be reindexed to the dimension of the main table.}
#' \item{\code{"response"}}{predictions after applying the inverse link function. If gp$lik.reindex2main = TRUE, 
#'						the predictions will be reindexed to the dimension of the main table.}
#' }
#'
#' @param !!!!!! upravit po nove definici - mozna udelat default NULL? CI numeric between 0 and 1, or \code{NULL}; confidence level for the confidence intervals to be calculated. If \code{NULL}, no standard errors
#' and confidence intervals will be calculated, which makes the prediction much faster (but then also no mean on the response scale will be
#' calculated - because for this, the standard error of the latent GP `f` is needed.)). Default is 0.95.
#' @param cov.fit logical; if \code{TRUE}, the covariance matrix of the predictions will also be returned (can be very memory consuming for large datasets!).
#' Default is \code{FALSE}.
#' @param link character, optional. In case of \code{type = "response"}, what should be the link function to take inverse of for calculating the derived quantity (response scale)? Character as passed to \code{make.link()}.
#' @param parname character, optional. In case of \code{type = "response"}, what should be the name of the parameter which is the derived quantity on the response scale?
#' @param maxn maximum maximum dataset size to fit at once (!!! the unit ???!!!); if the dataset is larger, it will be split into chunks of size \code{maxn}
#' before predictions are calculated. This parameter doesn't affect the result. It is only to speed-up the computation and use less memory by splitting the dataset
#' - the K(newdata,newdata) matrix is then smaller. So it only speeds it up for CI != NULL (in the backsolve function) for CI = NULL its only slowdown.
#'
#' @param pred.sims (currently unused) number of simulations to be used for calculating the response scale predictions and their CIs. Default is 100000.
#' Higher values give more accurate results, but are slower. Ignored if \code{CI = NULL}.
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
#	=> now is obsolete, we are using numerical integration. Whether or not it's used depends on CI.

# Note: weighting is done before predicting. Which is same as after, as it is all linear, but probably faster.

#
# groupMeans - 

# ... - passed to model_expand_predictions()



predict.gp <- function(gp, newdata = NULL, type = c('latent', 'terms', 'response'),	se.fit = FALSE, cov.fit = FALSE, CI = 0.95, 
						components = NULL, comp_missing = c("avg", "none"), w = NULL, groupMeans = NULL, 
						hyperpar = gpHyperparList(gp), 
						maxn = NULL, pred.sims = 100000, 
						Kx.cache = NULL, Kxx.cache = NULL, ...)  
{
	if (is.null(gp$fit))
		stop("Model object has not been fit yet: you need to call gpFit() first")
	comp_missing <- match.arg(comp_missing)
	need <- function (object, x) if (is.null(object[[x]])) stop("Model object is missing the `", x, "` element - try to call gpUnpack() on it")
	
	need(gp, "data")
	need(gp, "obsdata") #!!! na co?
	
	components <- validate_components(gp, components) # model components (components of covariance matrix)
	all_components_used <- all(names(gp$covComp) %in% components)
	
	type <- match.arg(type)
	
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
				if (!all(names(newdata[[tbl_name]]) %in% names(gp$obsdata[[tbl_name]]))) {
					stop("newdata$", tbl_name, " contains columns not present in training data")
				}
				if (!all(names(gp$obsdata[[tbl_name]] %in% names(newdata[[tbl_name]])))) {
					stop("newdata$", tbl_name, " is missing columns present in training data")
				}
            }
        }

		x_new <- gpDataPrepare(gp, newdata) # it will scale the data in exactly the way that training data set was		
	}
	
	# check CI
	if (!is.null(CI)) {
		if (CI >= 1 | CI <= 0) {
			stop("CI must be a number between 0 and 1, or NULL")
		}
		err <- qnorm( 1 - (1 - CI) / 2 )
		#se.fit <- TRUE
	}

	if (is.null(maxn)) maxn <- ceiling(gpDataSize(x_new, gp$GP_factor)  / 10) # we are splitting the newdata, not the training dataset; ceiling needed instead of round here
		
	# first, put up together the latent prediction
	mstart(id = "whole predict", mem_precise = TRUE)
	mstart(id = "whole pred", mem_precise = TRUE)
	pred <- pred(gp, x_new, same = same, hyperpar = hyperpar, components = components, comp_missing = comp_missing, 
		w = w, groupMeans = groupMeans, se.fit = se.fit, cov.fit = cov.fit, maxn = maxn, Kx.cache = Kx.cache, Kxx.cache = Kxx.cache
	)
	cat("pred() took ")
	mstop(id = "whole pred")
	if (cov.fit) {
		pred.cov <- pred$cov
		pred <- pred$pred
	}
	# pred is now always a matrix, with cols "f" and optionally also "f_SE" 
	
	if (type == "terms" || type == "response") {
		if (gp$lik.reindex2main && gp$GP_factor != "1") { 
			stopifnot(gpDataHasMainTable(x_new))	
			# reindex from the GP_factor to main table, as requested by the template
			fact_idx <- paste0(gp$GP_factor, "_idx")
			pred <- pred[x_new[[1]][[fact_idx]]]
			if (cov.fit)
				stop("type == terms or response used with cov.fit = TRUE and with lik.reindex2main = TRUE; would need to reindex the covariance matrix here, but sounds crazy, reconsider this case")
		}
		par <- c(hyperpar[[".lik"]], list(f = pred[,"f"]))
	}
	
	if (se.fit && !is.null(CI)) { # add confidence intervals
		upper <- pred[, "f"] + err * pred[, "f_SE"]
		lower <- pred[, "f"] - err * pred[, "f_SE"]
		dCI <- cbind(lower, upper)
		colnames(dCI) <- c(paste0("f_lower_", round(100 * CI), "CI"), #paste("lower ", round(100 * CI), "% CI", sep = ""),
						   paste0("f_upper_", round(100 * CI), "CI")) #paste("upper ", round(100 * CI), "% CI", sep = "")
		pred <- cbind(pred, dCI)							   
	}
	
	if (type == "terms") {
		terms <- gp$lik$terms(data = x_new, par) # we are calling $terms and not $stage(stage = 1), because we want just the terms separately
		if (!is.null(terms))
			pred <- cbind(pred, as.matrix(terms)) # we want to keep the return value as matrix, as it has always been
	}
	else if (type == "response") {
		pred <- as.matrix(gp$lik$stages(data = x_new, par, stages = 1:2)) # in this case, we don't return f and f_SE and the other stuff at the moment... we could though... to reconsider
			# we want to keep the return value as matrix, as it has always been
	}
	cat("returning memory - gc() took ")
	mstart(id = "gc")
	gc() # Return the memory. Important! :-) 
	mstop(id = "gc")
	cat("whole predict() took ")
	mstop(id = "whole predict")
	if (!cov.fit)
		pred
	else
		list(pred = pred, cov = pred.cov)
}