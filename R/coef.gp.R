#' Extract fitted Gaussian-process hyperparameters
#'
#' Returns the hyperparameter values used for a fitted \code{gp} model together
#' with the marginal-log-likelihood, prior, and posterior gradients.
#'
#' The column \code{fitted value (h)} is on the hyperparameter scale used by the
#' covariance or likelihood function. It is not on the optimizer scale.
#'
#' The column \code{gradient: d log(marginal likelihood) / d h} is
#' \eqn{\partial \log p(y \mid h) / \partial h}, evaluated at the displayed
#' fitted value. More precisely, it is the gradient of the approximate marginal
#' log likelihood stored in \code{object$fit$h_grads}.
#'
#' The column \code{gradient: d log(prior) / d h} is
#' \eqn{\partial \log p(h) / \partial h}, evaluated from the hyperparameter
#' priors configured in \code{object$hyperpar}. It is displayed whether or not
#' the prior was used when fitting the model.
#'
#' The column \code{gradient: d log(posterior) / d h} is
#' \eqn{\partial \log p(h \mid y) / \partial h}, calculated as the sum of the
#' marginal-log-likelihood and log-prior gradients. Thus:
#' \itemize{
#'	\item all three gradients are with respect to the displayed
#'		hyperparameter \eqn{h}, not the transformed optimizer parameter;
#'	\item the marginal-log-likelihood gradient has the opposite sign from the
#'		gradient of \code{object$fit$nlml};
#'	\item the posterior gradient has the opposite sign from the gradient of
#'		the negative log posterior minimized when the prior is used for
#'		optimization.
#' }
#' In each gradient column, a positive value means that locally increasing
#' \eqn{h} increases the corresponding log quantity; a negative value means
#' that locally decreasing \eqn{h} increases it.
#'
#' Fixed hyperparameters are included to show the complete fitted
#' hyperparameter set, but their gradients are \code{NA}, because
#' \code{object$fit$h_grads} contains entries only for non-fixed
#' hyperparameters. The marginal-log-likelihood and posterior gradients can
#' also be \code{NA} when marginal-log-likelihood gradient computation was
#' disabled.
#'
#' @param object A fitted object of class \code{gp}.
#' @param ... Additional arguments, currently ignored.
#'
#' @return A data frame with one row per hyperparameter, in the same order as
#' \code{object$hyperpar}.
#'
#' @method coef gp
#' @export
coef.gp <- function(object, ...)
{
	fit <- object[["fit"]]
	if (is.null(fit))
		stop("Model object has not been fit yet: you need to call gpFit() first")

	h <- fit[["h"]]
	h_grads <- fit[["h_grads"]]
	if (is.null(h))
		stop("The fitted model does not contain hyperparameter values in `object$fit$h`.")
	if (is.null(h_grads))
		stop("The fitted model does not contain hyperparameter gradients in `object$fit$h_grads`.")

	nonfixed <- !object$hyperpar$fixed
	if (length(h) != sum(nonfixed))
		stop("The number of fitted hyperparameters in `object$fit$h` does not match the non-fixed rows of `object$hyperpar`.")
	if (length(h_grads) != sum(nonfixed))
		stop("The number of gradients in `object$fit$h_grads` does not match the non-fixed rows of `object$hyperpar`.")

	fitted_value <- object$hyperpar$value
	fitted_value[nonfixed] <- h
	marginal_likelihood_gradient <- rep(NA_real_, nrow(object$hyperpar))
	prior_gradient <- rep(NA_real_, nrow(object$hyperpar))
	posterior_gradient <- rep(NA_real_, nrow(object$hyperpar))

	marginal_likelihood_gradient[nonfixed] <- h_grads
	prior_gradient[nonfixed] <- gpPriorGradient(object, h)
	posterior_gradient[nonfixed] <-
		marginal_likelihood_gradient[nonfixed] + prior_gradient[nonfixed]

	data.frame(
		component = object$hyperpar$component,
		hyperpar = object$hyperpar$hyperpar,
		i = object$hyperpar$i,
		table = object$hyperpar$table,
		var = object$hyperpar$var,
		fixed = object$hyperpar$fixed,
		"fitted value (h)" = fitted_value,
		"gradient: d log(marginal likelihood) / d h" = marginal_likelihood_gradient,
		"gradient: d log(prior) / d h" = prior_gradient,
		"gradient: d log(posterior) / d h" = posterior_gradient,
		check.names = FALSE
	)
}
