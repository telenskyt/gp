#' Convert latent GP predictions to the requested prediction type
#'
#' `gpPredictFromLatent()` is an internal helper used by prediction functions
#' that have already computed latent predictions for `f`. It applies the common
#' `type = "latent"`, `"terms"`, and `"response"` prediction logic used by
#' `predict.gp()`.
#'
#' @param gp Object of class `gp`.
#' @param pred Matrix or data frame with column `f`, and optionally `f_SE`.
#'	The rows must initially be in the dimension of `gp$GP_factor`.
#' @param data Object of class `gpData` to pass to the likelihood template.
#'	The caller must provide data on the scale expected by the template.
#' @param type Character; prediction type. Same meaning as in `predict.gp()`.
#' @param hyperpar List of hyperparameter values, usually from `gpHyperparList(gp)`.
#'	When likelihood	hyperparameters are present in the model (fixed or fitted), 
#'  it must already contain the `.lik` component with correct hyperparameter values.
#' @param se.fit Logical; should standard errors be present in the returned
#'	prediction matrix? If `TRUE` and `pred` has no `f_SE` column, `pred.cov`
#'	must be supplied.
#' @param cov.fit Logical; should the latent covariance matrix be returned?
#' @param pred.cov Optional covariance matrix of latent predictions, in the same
#'	row order as `pred` before any likelihood-template reindexing.
#' @param conf.int Logical; should latent-scale confidence interval columns be
#'	added? Same meaning as in `predict.gp()`.
#' @param conf.level Numeric confidence level used when `conf.int = TRUE`.
#'
#' @return If `cov.fit = FALSE`, a prediction matrix. If `cov.fit = TRUE`, a list
#'	with components `pred` and `cov`, matching `predict.gp()`.
#'
#' For `type = "latent"`, predictions remain in the dimension of `gp$GP_factor`.
#' For `type = "terms"` or `type = "response"`, predictions are reindexed to the
#' main table when `gp$lik.reindex2main` requires it. In that reindexed case,
#' `cov.fit = TRUE` is currently not supported, matching `predict.gp()`.
#'
#' @keywords internal
gpPredictFromLatent <- function(gp, pred, data, type = c("latent", "terms", "response"),
	hyperpar = gpHyperparList(gp), se.fit = FALSE, cov.fit = FALSE,
	pred.cov = NULL, conf.int = FALSE, conf.level = 0.95)
{
	type <- match.arg(type)
	pred <- as.matrix(pred)

	if (is.null(colnames(pred)) || !"f" %in% colnames(pred))
		stop("pred must contain a column named 'f'")
	if (!is.logical(se.fit) || length(se.fit) != 1L || is.na(se.fit))
		stop("se.fit must be TRUE or FALSE")
	if (!is.logical(cov.fit) || length(cov.fit) != 1L || is.na(cov.fit))
		stop("cov.fit must be TRUE or FALSE")
	if (!is.logical(conf.int) || length(conf.int) != 1L || is.na(conf.int))
		stop("conf.int must be TRUE or FALSE")
	if (!is.numeric(conf.level) || length(conf.level) != 1L || is.na(conf.level) || conf.level >= 1 || conf.level <= 0)
		stop("conf.level must be a number between 0 and 1")

	if (cov.fit && is.null(pred.cov))
		stop("pred.cov must be supplied when cov.fit = TRUE")
	if ((se.fit || cov.fit || conf.int) && !"f_SE" %in% colnames(pred)) {
		if (is.null(pred.cov))
			stop("pred must contain column 'f_SE', or pred.cov must be supplied")
		pred <- cbind(pred, f_SE = sqrt(diag(pred.cov)))
	}
	# pred is now always a matrix, with cols "f" and optionally also "f_SE"

	if (type == "terms" || type == "response") {
		# do the same as in gpGetParForLikTemplate(), but reindex not only f, but also pred!
		if (gp$lik.reindex2main && gp$GP_factor != "1") {
			stopifnot(gpDataHasMainTable(data))
			# reindex from the GP_factor to main table, as requested by the template
			fact_idx <- paste0(gp$GP_factor, "_idx")
			pred <- pred[data[[1]][[fact_idx]],,drop=FALSE]
			if (cov.fit)
				stop("type == terms or response used with cov.fit = TRUE and with lik.reindex2main = TRUE; would need to reindex the covariance matrix here, but sounds crazy, reconsider this case")
		}
		par <- c(hyperpar[[".lik"]], list(f = pred[,"f"]))
	}

	if (conf.int) { # add confidence intervals
		err <- qnorm( 1 - (1 - conf.level) / 2 )
		upper <- pred[, "f"] + err * pred[, "f_SE"]
		lower <- pred[, "f"] - err * pred[, "f_SE"]
		dCI <- cbind(lower, upper)
		colnames(dCI) <- c(paste0("f_lower_", round(100 * conf.level), "CI"),
						   paste0("f_upper_", round(100 * conf.level), "CI"))
		pred <- cbind(pred, dCI)
	}

	if (type == "terms") {
		terms <- gp$lik$terms(data = data, par) # we are calling $terms and not $stage(stage = 1), because we want just the terms separately
		if (!is.null(terms))
			pred <- as.matrix(bind_cols(pred, as.matrix(terms))) # we want to keep the return value as matrix, as it has always been
	}
	else if (type == "response") {
		pred <- as.matrix(gp$lik$stages(data = data, par, stages = 1:2)) # in this case, we don't return f and f_SE and the other stuff at the moment... we could though... to reconsider
			# we want to keep the return value as matrix, as it has always been
	}

	stopifnot(is.matrix(pred))
	if (!cov.fit)
		pred
	else
		list(pred = pred, cov = pred.cov)
}
