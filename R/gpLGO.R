
#' Calculate the approximation of the LGO-CV (leave-group-out cross-validation) and evaluate it.
#'
#' @param gp GP model object
#' @param fold.col either a name of a column in the main table, or a vector along factor \code{fold.fact}.
#'		The column (or the supplied vector) must be a vector of integers from \code{1} to \code{N} (\code{N} being the number of folds),
#'		specifying the number of the cross-validation fold the given record belongs to.
#' @param fold.fact a factor along which the cross-validation folds (\code{fold.col}) are specified. Factor \code{"1"} means the folds are specified for the rows of the main table.
#'		Note that if the GP dimension is given by some real grouping factor (i.e. \code{gp$GP_factor != "1"}), then the \code{fold.fact} must be that factor.
#' @param pred.null optional predictions from the null model (\code{predict.lmFit(type = "terms")}). If these are available from the actual exact lmFitCV call,
#'		it is better to provide them this way. Otherwise the null model predictions will be calculated just by fitting all data as training, which will result in a slightly
#'		better null model (and a slightly lower pseudo-R2 statistics then).
#' @param lmFit.options list of options to be passed to \code{\link{lmFit}} when fitting the intercept-only null model. Only used when \code{pred.null = NULL}.
#' @param lmFit.pred.options list of options to be passed to the \code{predict.lmFit} method when predicting from the null model. Only used when \code{pred.null = NULL}.
#' @param hyperpar (optional) lists of hyperparameter values, as returned by \code{gpHyperparList()}; default is to use the hyperparameters from the \code{gp} object.
#'
#' @returns the gp model object with approximate LGO-CV predictions and statistics in \code{gp$fit$LGOCV}.
gpLGO <- function(gp, fold.col, fold.fact = "1",
	pred.null = NULL,
	lmFit.options = list(lik.hyperpar = "fix"),
	lmFit.pred.options = list(type = "terms"),
	hyperpar = gpHyperparList(gp))
{
	if (!is.null(lmFit.options)) {
		if (any(c("gp", "formula", "data") %in% names(lmFit.options)))
			stop("lmFit.options must not contain gp, formula, or data argument")
	}

	gp$fit$LGOCV <- gpLGOpred(gp, fold.col = fold.col, fold.fact = fold.fact)

	mstart(TRUE, id = "LGO_pred_dens")
	LGOev <- pred_dens_logit__int(gp, pred = gp$fit$LGOCV)
	cat("predictive density calculation took ")
	mstop(TRUE, id = "LGO_pred_dens")

	pred <- as.matrix(data.frame(f = as.numeric(gp$fit$LGOCV$f)))

	pred <- cbind(pred, f_SE = sqrt(diag(gp$fit$LGOCV$cov)))

	pred <- gpPredictFromLatent(gp = gp, pred = pred, data = gp$data, type = "terms",
		hyperpar = hyperpar, se.fit = TRUE,
		pred.cov = gp$fit$LGOCV$cov)

	if (is.null(pred.null)) { # get null model predictions
		nullFit.args <- lmFit.options
		nullFit.args$gp <- NULL
		nullFit.args$formula <- NULL
		nullFit.args$data <- NULL
		nullm <- do.call(function(...) lmFit(gp, formula = ~1, data = gp$data, ...), nullFit.args)
		pred.null <- do.call(function(...) predict(nullm, newdata = gp$data, ...), lmFit.pred.options)
	}

	gp$fit$LGOCV$stats <- gp$lik$predPerf(data = gp$data, pred = as.data.frame(pred), pred.null = as.data.frame(pred.null))
	gp$fit$LGOCV$stats$NLPD <- LGOev$value
	gp
}
