# cov.NN(x1, x2, l)
# Sum of Neural Network Covariance Functions, a la Kallasvuo, Vanhatalo & Veneranta 2017, dx.doi.org/10.1139/cjfas-2016-0008
#
# x1 - matrix of covariates for the training dataset (rows - data records, cols - covariates)
# x2 - matrix of covariates for the predicted dataset (rows - data records, cols - covariates)
#
# argument variants:
# 1) no x2
# 3) with x2

# output: 
#	- covariance matrix of dimensions nrow(x1) x nrow(x1) (variants 1)
#	- covariance matrix of dimensions nrow(x1) x nrow(x2) (variants 3)

# Internal version cov.NN.additive.naive.b_correct_arguments.R

#source("cov.NN.R")

cov.NN.add <- function(x1, x2 = NULL, sigma2_int, sigma2_slope) 
{
	stopifnot(length(sigma2_int) == ncol(x1))
	stopifnot(length(sigma2_slope) == ncol(x1))
	stopifnot(is.null(x2) || ncol(x1) == ncol(x2))
	K <- matrix(0, nrow = nrow(x1), ncol = if (is.null(x2)) nrow(x1) else nrow(x2))
	for (i in 1:ncol(x1)) {
		if (is.null(x2)) {
			K <- K + cov.NN(x1[,i,drop = FALSE], sigma2_diag = c(sigma2_int[i], sigma2_slope[i]))
		} else {
			K <- K + cov.NN(x1[,i,drop = FALSE], x2[,i,drop = FALSE], sigma2_diag = c(sigma2_int[i], sigma2_slope[i]))
		}
		gc() # this is needed to keep memory usage reasonably low, in parallel processing of big matrices...
	}
	K
}
