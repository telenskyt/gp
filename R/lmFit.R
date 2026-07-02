# The  version generated with OpenAI Codex, after a long conversation. I instructed Codex to refactor it from the code 
# specific for the Lapwing case study (nestsurv-lm.R). Now I regret a bit, the code is more complicated than I would do it,
# even after many requests for simplification...
# !!! Perhaps the lik.hyperpar = "optimize" branch might be buggy.
# 
#
#' Fit a linear model on the GP model object
#'
#' Fits a linear model, where the latent Gaussian process \code{f}
#' is replaced by a linear predictor constructed from \code{formula}.
#'
#' @param gp gp model object; its likelihood template will be used
#' @param formula the linear model formula to replace the latent Gaussian process \code{f}
#' @param data object of class \code{gpData}; data used for the likelihood and the linear predictor. See the details below.
#' @param table character - name of the table in \code{data} that is going to be used to construct the model matrix for 
#'	the linear model. This table must be of the dimension of \code{gp$GP_factor}.
#' @param beta.name character; name of the linear predictor coefficient parameter
#' @param lik.hyperpar character; if there are likelihood hyperparameters, should they be fixed or optimized?
#' \describe{
#'	\item{\code{"fix"}}{the likelihood hyperparameters will be fixed at value specified by \code{lik.fix}}
#'	\item{\code{"optimize"}}{the likelihood hyperparameters will be optimized along with the parameters in \code{formula}}
#' }
#' @param lik.fix when \code{lik.hyperpar = "fix"}, to what value should the likelihood hyperparameters be fixed?
#' \describe{
#'	\item{single numeric value}{all likelihood hyperparameters will be set to this value}
#'	\item{\code{"value"}}{the likelihood hyperparameters will be set to values in the column \code{"value"} in the hyperparameter table (\code{gp$hyperpar$value})}
#' }
#' @param lik.prefix character; prefix added to likelihood hyperparameter names in \code{summary()}
#' @param silent logical; passed to \code{RTMB::MakeADFun}
#' @param optCtrl list; passed to \code{nlminb(control = )}
#'
#' @details Note that the model is linear in the \code{formula} specified. If \code{lik.hyperpar = "optimize"}, 
#' the model is not necessarily linear in these hyperparameters; this is given by the likelihood template defined in the \code{gp} object.
#'
#' The table used for the linear model must have the same dimension as \code{gp$GP_factor}, because the linear predictor replaces
#' the latent Gaussian process \code{f} directly. If \code{gp$GP_factor == "1"}, the table must have no grouping factor.
#'
#' Note that no scaling (covariate standardization) is done neither in this function nor in \code{predict}. It is up to the user 
#' to provide the correct \code{data} argument with the desired scaling.
#'
#' This function is used in the cross-validation process to fit the null model for the baseline comparison.
#'
#' @return Object of class \code{lmFit}
#' @export
lmFit <- function(gp, formula, data, table, beta.name = "beta_lm", lik.hyperpar = c("fix", "optimize"), lik.fix = NULL, lik.prefix = ".lik.", silent = TRUE, optCtrl = list(iter.max=300, eval.max=400))
{
	# Check basic inputs.
	stopifnot(class(gp) == "gp")
	stopifnot(is(gp$lik, "likTempPhased"))
	stopifnot(class(data) == "gpData")
	stopifnot(is.character(table) && length(table) == 1)
	stopifnot(is.character(beta.name) && length(beta.name) == 1 && nchar(beta.name) > 0)
	lik.hyperpar <- match.arg(lik.hyperpar)
	stopifnot(is.character(lik.prefix) && length(lik.prefix) == 1)
	lmFitCheckTableFactor(data, table, gp, caller = "lmFit()", data.name = "data")

	# Build the model matrix for the linear predictor.
	data_table <- as.data.frame(data[[table]])
	x <- model.matrix(formula, data_table)
	if (nrow(x) != nrow(data[[table]]))
		stop("model.matrix() changed the number of rows from ", nrow(data[[table]]), " to ", nrow(x), ".")

	# Work on a copy of the GP object without fitted state.
	gp_lm <- gp
	gp_lm$fit <- NULL

	# Decide which likelihood hyperparameters are optimized and which fixed values are used.
	lik.ind <- gp_lm$hyperpar$component == ".lik"
	# check and process the lik.fix argument if it applies
	if (!any(lik.ind) && lik.hyperpar == "optimize")
		stop("there are no likelihood hyperparameters to optimize")
	if (lik.hyperpar == "fix" && any(lik.ind)) { # there are likelihood hyperparameters
		if (is.null(lik.fix))
			stop("the likelihood template has hyperparameters, and you specified lik.hyperpar = 'fix'. Please specify how should they be fixed with the lik.fix argument.")
		else if (is.character(lik.fix) && lik.fix == "value") {
			# everything is fine
		} else if (is.numeric(lik.fix) && length(lik.fix) == 1) {
			gp_lm$hyperpar$value[lik.ind] <- rep(lik.fix, length.out = sum(lik.ind))
		} else {
			stop("numeric lik.fix needs to have length exactly 1 (a single number)")
		}
	}
		
	gp_lm$hyperpar$fixed <- TRUE
	if (lik.hyperpar == "optimize")
		gp_lm$hyperpar$fixed[lik.ind] <- FALSE
	lik.optim.ind <- !gp_lm$hyperpar$fixed

	# Warn about a duplicate intercept in the linear predictor and optimized likelihood parameters.
	if (lik.hyperpar == "optimize" && "(Intercept)" %in% colnames(x) &&
		any(lik.optim.ind & !is.na(gp_lm$hyperpar$var) & gp_lm$hyperpar$var == "(Intercept)"))
		warning(
			"The linear model formula and optimized likelihood hyperparameters both contain an intercept; ",
			"this can make the model non-identifiable."
		)

	# Build the parameter list passed to RTMB.
	par <- structure(list(structure(numeric(ncol(x)), names = colnames(x))), names = beta.name)
	lik.par.names <- character(0)
	if (lik.hyperpar == "optimize") {
		lik.par <- gpHyperparList(gp_lm, "start")[[".lik"]]
		lik.par.names <- names(lik.par)
		if (length(lik.par) > 0)
			par <- c(par, lik.par)
	}
	if (anyDuplicated(names(par)))
		stop("Parameter names are not unique.")

	# Build the RTMB objective function.
	cmb <- function(f, gp, data, x, beta.name, lik.par.names) function(p) f(p, gp, data, x, beta.name, lik.par.names)

	obj <- MakeADFun(cmb(lmFitNll, gp_lm, data, x, beta.name, lik.par.names), par, silent = silent)

	# Build names for coef()/summary() from the exact parameter list passed to
	# RTMB, rather than reconstructing the order from the hyperparameter table.
	# This keeps the displayed names tied to the optimized vector order.
	coef.names <- names(unlist(par, use.names = TRUE))
	beta.prefix <- paste0(beta.name, ".")
	beta.ind <- startsWith(coef.names, beta.prefix)
	coef.names[beta.ind] <- substring(coef.names[beta.ind], nchar(beta.prefix) + 1)
	coef.names[!beta.ind] <- paste0(lik.prefix, coef.names[!beta.ind])

	# RTMB flattens list parameters by repeating the top-level list name for
	# each scalar entry, e.g. names(list(beta_lm = c(a = 0, b = 0))) becomes
	# c("beta_lm", "beta_lm"). Check that this is exactly the same flattening
	# as our input list, because coef.names are built from that list order.
	if (!identical(names(obj$par), rep(names(par), lengths(par))) || length(coef.names) != length(obj$par))
		stop("RTMB parameter order does not match the parameter list.")

	# Optimize the objective.
	fit <- nlminb(obj$par, obj$fn, obj$gr, control = optCtrl)
	if (!is.null(fit$convergence) && fit$convergence != 0) {
		stop("Model convergence problem; ",
			fit$message, ". Consider increasing iter.max and eval.max in the optCtrl argument.")
	}

	lik.par.fit <- NULL
	if (lik.hyperpar == "optimize" && any(lik.optim.ind))
		lik.par.fit <- obj$env$parList(fit$par)[lik.par.names]

	# Compute the standard errors and p-values.
	if (!silent)
		cat("sdreport:\n")
	sdr <- sdreport(obj)

	# Return the fitted object.
	object <- list(
		call = match.call(),
		formula = formula,
		table = table,
		data = data,
		x = x,
		gp = gp_lm,
		obj = obj,
		fit = fit,
		sdr = sdr,
		beta.name = beta.name,
		lik.hyperpar = lik.hyperpar,
		lik.par = lik.par.fit,
		lik.prefix = lik.prefix,
		coef.names = coef.names
	)
	class(object) <- "lmFit"
	object
}

# Evaluate the negative log likelihood for the linearized model.
lmFitNll <- function(par, gp, data, x, beta.name, lik.par.names)
{
	# Compute the linear predictor from the model matrix and linear coefficients.
	f <- (x %*% par[[beta.name]])[,1]

	hyperpar <- gpHyperparList(gp)
	if (length(lik.par.names) > 0)
		hyperpar[[".lik"]] <- par[lik.par.names]

	# Let the standard likelihood helper handle .lik and possible GP-factor-to-main-table reindexing.
	lik.par <- gpGetParForLikTemplate(gp, f, data, hyperpar)
	gp$lik$nll(data, lik.par)
}

lmFitCheckTableFactor <- function(data, table, gp, caller, data.name)
{
	if (!table %in% names(data))
		stop("Table `", table, "` is missing in `", data.name, "`.")
	table.fact <- attr(data[[table]], "fact")
	if (is.null(table.fact))
		table.fact <- "1"
	if (table.fact != gp$GP_factor)
		stop(
			caller, " requires table `", table, "` in `", data.name, "` to be of the dimension of `gp$GP_factor`. ",
			"Table `", table, "` has factor `", table.fact, "`, but `gp$GP_factor` is `", gp$GP_factor, "`."
		)
	invisible(table.fact)
}

#' @method terms lmFit
#' @export
terms.lmFit <- function(object, ...)
{
	terms(object$formula)
}

#' @method coef lmFit
#' @export
coef.lmFit <- function(object, ...)
{
	sum <- summary(object$sdr, p.value = TRUE)
	rownames(sum) <- object$coef.names
	if (!all(is.finite(sum))) {
		warning("summary.sdreport() produced non-numeric values. One of the reasons might be a singular matrix.", 
			if (object$lik.hyperpar == "optimize") 
				" This can e.g. happen if there are two intercepts, one in the linear predictor and the other in the (now optimized) likelihood hyperparameters."
			else
				""
		)
	}
	sum
}

#' @method predict lmFit
#' @export
predict.lmFit <- function(object, newdata = NULL, se.fit = FALSE, ...)
{
	if (se.fit == TRUE && object[["lik.hyperpar"]] == "optimize")
		stop("se.fit = TRUE for lik.hyperpar = 'optimize' not implemented at the moment")
	# Use the training table unless the caller supplied prediction data.
	if (is.null(newdata)) {
		newdata <- object$data
		data.name <- "object$data"
	} else {
		stopifnot(inherits(newdata, "gpData"))
		data.name <- "newdata"
	}
	lmFitCheckTableFactor(newdata, object$table, object$gp, caller = "predict.lmFit()", data.name = data.name)
	data_table <- as.data.frame(newdata[[object$table]])

	# Build the prediction model matrix.
	A <- model.matrix(object$formula, data_table)
	if (nrow(A) != nrow(data_table))
		stop("model.matrix() changed the number of rows from ", nrow(data_table), " to ", nrow(A), ".")
	if (!identical(colnames(A), colnames(object$x)))
		stop("Prediction model matrix columns do not match the fitted model matrix.")

	# Extract fitted RTMB parameters and compute the linear predictor.
	par <- object$obj$env$parList(object$fit$par)
	f <- (A %*% par[[object$beta.name]])[,1]

	# Return only the linear predictor unless standard errors were requested.
	if (!se.fit)
		return(f)

	# Compute standard errors for the linear predictor from the beta covariance block.
	sigma <- solve(object$obj$he(object$fit$par))
	beta.ind <- seq_len(ncol(object$x))
	sigma.beta <- sigma[beta.ind, beta.ind, drop = FALSE]
	f.SE <- sqrt(rowSums((A %*% sigma.beta) * A))
	data.frame(f = f, f_SE = f.SE)
}

# Count the parameters optimized by RTMB.
lmFitDf <- function(object)
{
	length(object$fit$par)
}

#' @method nobs lmFit
#' @export
nobs.lmFit <- function(object, use.fallback = FALSE, ...)
{
	# Prefer the number of likelihood observations reported by the RTMB object.
	rep <- object$obj$report(par = object$fit$par)
	if (!is.null(rep$y))
		return(length(rep$y))

	# Fall back to model-matrix rows only when requested by the caller.
	if (use.fallback)
		return(nrow(object$x))

	stop("Cannot determine number of observations: `obj$report()` did not return `y`.", call. = FALSE)
}

#' @method logLik lmFit
#' @export
logLik.lmFit <- function(object, ...)
{
	structure(
		-object$fit$objective,
		df = lmFitDf(object),
		nobs = stats::nobs(object, use.fallback = TRUE),
		class = "logLik"
	)
}

#' @method AIC lmFit
#' @export
AIC.lmFit <- function(object, ..., k = 2)
{
	2*object$fit$objective + k*lmFitDf(object)
}

#' @method extractAIC lmFit
#' @export
extractAIC.lmFit <- function(object, scale, k = 2, ...)
{
	c(lmFitDf(object), AIC.lmFit(object, k = k))
}

#' @method summary lmFit
#' @export
summary.lmFit <- function(object, ...)
{
	cat("formula: ")
	print(object$formula)
	cat("table: ", object$table, "\n\n", sep = "")

	printCoefmat(coef(object), ...)

	cat("\nnegative log likelihood: ", object$fit$objective, "\n", sep = "")
	cat("AIC: ", stats::AIC(object), "\n", sep = "")
	invisible(object)
}

#' @method print lmFit
#' @export
print.lmFit <- function(x, ...)
{
	summary.lmFit(x, ...)
}
