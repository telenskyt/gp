# cov.NN(x1, x2, l)
# Sum of Neural Network Covariance Functions, a la Kallasvuo, Vanhatalo & Veneranta 2017, dx.doi.org/10.1139/cjfas-2016-0008
#
# x1 - matrix of covariates for the training dataset (rows - data records, cols - covariates)
# x2 - matrix of covariates for the predicted dataset (rows - data records, cols - covariates)
#
# argument variants:
# 1) no x2
# 3) with x2 (not now, currently not used)

# output: 
#	- covariance matrix of dimensions nrow(x1) x nrow(x1) (variants 1)
#	- covariance matrix of dimensions nrow(x1) x nrow(x2) (variants 3)

# Based on the fastest version, cov.NN.additive.naive.b_correct_arguments.R

#source("cov.NN.d1.R")

# j = 1 intercept, j = 2 slope
cov.NN.add.d1 <- function(x1, sigma2_int, sigma2_slope, der_wrt, der_wrt_i) 
{
	stopifnot(length(sigma2_int) == ncol(x1))
	stopifnot(length(sigma2_slope) == ncol(x1))
	stopifnot(der_wrt_i <= length(sigma2_int) && der_wrt_i >= 1)
	if (der_wrt == "sigma2_int")
		i <- 1
	else if (der_wrt == "sigma2_slope")
		i <- 2
	else
		stop("Invalid value of der_wrt ('", der_wrt, "')")

	cov.NN.d1(x1[,der_wrt_i,drop = FALSE], x2 = NULL, sigma2_diag = c(sigma2_int[der_wrt_i], sigma2_slope[der_wrt_i]), der_wrt = "sigma2_diag", der_wrt_i = i)
}
