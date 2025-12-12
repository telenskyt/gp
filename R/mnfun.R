
# hyperpar - list of hyperparameters
#' Compute mean of the Gaussian Process
#' 
#' Compute mean of the Gaussian Process along given dataset. Currently just returns vector of zeros.
#'
#' @param gp object of class \code{gp}.
#' @param data data along which to generate the gaussian process mean function. 
#' @param hyperpar (optional) lists of hyperparameter values, as returned by \code{gpHyperparList()}; default is to use the hyperparameters from the \code{gp} object.
#' 
#' @return Vector of the length of \code{gpDataSize(data, gp$GP_factor)}.

mnfun <- function (gp, data = gp$data, hyperpar = NULL)
{
	rep(0, gpDataSize(data, gp$GP_factor))
}
