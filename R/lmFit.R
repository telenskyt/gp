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
#' @param formula the linear model formula
#' @param data object of class \code{gpData}; data used for the likelihood and the linear predictor. See the details below.
#' @param table character - name of the table in \code{data} that is going to be used to construct the model matrix for 
#'	the linear model
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
	if (!table %in% names(data))
		stop("Table `", table, "` is missing in the data.")

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
	if (lik.hyperpar == "optimize") {
		lik.par <- gpHyperparExportVector(gp_lm, "start")
		if (length(lik.par) > 0)
			par <- c(par, setNames(list(lik.par), ".lik"))
	}
	if (anyDuplicated(names(par)))
		stop("Parameter names are not unique.")

	# Build coefficient names shown by coef() and summary().
	coef.names <- colnames(x)
	if (lik.hyperpar == "optimize" && any(lik.optim.ind)) {
		lik.names <- gp_lm$hyperpar$hyperpar[lik.optim.ind]
		lik.var <- gp_lm$hyperpar$var[lik.optim.ind]
		lik.names[!is.na(lik.var)] <- paste(lik.names[!is.na(lik.var)], lik.var[!is.na(lik.var)], sep = ".")
		coef.names <- c(coef.names, paste0(lik.prefix, lik.names))
	}

	# Build the RTMB objective function.
	cmb <- function(f, gp, data, x, table, beta.name) function(p) f(p, gp, data, x, table, beta.name)

	obj <- MakeADFun(cmb(lmFitNll, gp_lm, data, x, table, beta.name), par, silent = silent)

	# Optimize the objective.
	fit <- nlminb(obj$par, obj$fn, obj$gr, control = optCtrl)
	if (!is.null(fit$convergence) && fit$convergence != 0) {
		stop("Model convergence problem; ",
			fit$message, ". Consider increasing iter.max and eval.max in the optCtrl argument.")
	}

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
		lik.prefix = lik.prefix,
		coef.names = coef.names
	)
	class(object) <- "lmFit"
	object
}

# Evaluate the negative log likelihood for the linearized model.
lmFitNll <- function(par, gp, data, x, table, beta.name)
{
	# Compute the linear predictor from the model matrix and linear coefficients.
	f <- (x %*% par[[beta.name]])[,1]

	# Reindex the linear predictor to the dimension expected by the likelihood.
	f <- lmFitReindexF(gp, f, data, table)

	# Remove the linear coefficients; the likelihood only needs hyperparameters and f.
	par[[beta.name]] <- NULL

	# Select the likelihood hyperparameter rows.
	lik.ind <- gp$hyperpar$component == ".lik"
	lik.hy <- gp$hyperpar[lik.ind, , drop = FALSE]

	# Use RTMB values when .lik is optimized; otherwise use fixed values from the table.
	if (!is.null(par[[".lik"]]))
		lik.value <- par[[".lik"]]
	else
		lik.value <- gp$hyperpar$value[lik.ind]

	# Convert likelihood hyperparameters from table rows to the likelihood list.
	lik.par <- setNames(lapply(unique(lik.hy$hyperpar), function(hyperpar) {
		lik.value[lik.hy$hyperpar == hyperpar]
	}), unique(lik.hy$hyperpar))

	# Add the linear predictor and evaluate the likelihood.
	lik.par$f <- f
	gp$lik$nll(data, lik.par)
}

# Determine the factor dimension expected by the likelihood.
lmFitLikFactor <- function(gp)
{
	if (gp$lik.reindex2main && gp$GP_factor != "1")
		return("1")
	gp$GP_factor
}

# Determine the factor dimension of the table used by the linear predictor.
lmFitTableFactor <- function(data, table)
{
	fact <- attr(data[[table]], "fact")
	if (is.null(fact))
		fact <- "1"
	fact
}

# Reindex the linear predictor to the dimension expected by the likelihood.
lmFitReindexF <- function(gp, f, data, table)
{
	source_fact <- lmFitTableFactor(data, table)
	target_fact <- lmFitLikFactor(gp)

	if (source_fact == target_fact) {
		target_size <- gpDataSize(data, target_fact)
		if (length(f) != target_size)
			stop("The linear predictor has length ", length(f), ", but the likelihood expects length ", target_size, ".")
		return(f)
	}

	if (target_fact == "1" && source_fact != "1") {
		if (!gpDataHasMainTable(data))
			stop("Cannot reindex the linear predictor from factor `", source_fact, "` to the main table: the data has no main table.")
		fact_idx <- paste0(source_fact, "_idx")
		if (!fact_idx %in% colnames(data[[1]]))
			stop("Cannot reindex the linear predictor: column `", fact_idx, "` is missing in the main table.")
		return(f[data[[1]][[fact_idx]]])
	}

	stop("Cannot reindex the linear predictor from factor `", source_fact, "` to factor `", target_fact, "`.")
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
	# Use the original table unless the caller supplied prediction data.
	if (is.null(newdata)) {
		data_table <- as.data.frame(object$data[[object$table]])
	} else if (inherits(newdata, "gpData")) {
		if (!object$table %in% names(newdata))
			stop("Table `", object$table, "` is missing in `newdata`.")
		data_table <- as.data.frame(newdata[[object$table]])
	} else {
		data_table <- as.data.frame(newdata)
	}

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
