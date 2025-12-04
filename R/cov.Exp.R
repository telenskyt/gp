# cov.Exp(x1, x2, l) ; (~ Matern with v = 1/2)
#
# x1 - matrix of covariates for the training dataset (rows - data records, cols - covariates)
# x2 - matrix of covariates for the predicted dataset (rows - data records, cols - covariates)
#
# ls - length-scales = as defined by Rasmussen & Williams 2006, which are length-scales at the scale of the covariates
#		note that in the paper Golding & Purse 2016 is an opposite thing, the l is defined as a sqrt of the Rasmussen & Williams 2006 length-scale (super confusion there :))
#
# argument variants:
# 1) no x2
# 3) with x2

# output: 
#	- covariance matrix of dimensions nrow(x1) x nrow(x1) (variant 1)
#	- covariance matrix of dimensions nrow(x1) x nrow(x2) (variant 3)
#
# Reference: it is used here as a spatial covariance function: Kallasvuo, M., Vanhatalo, J., & Veneranta, L. (2016). Modeling the spatial distribution of larval fish abundance provides essential information for management. Canadian Journal of Fisheries and Aquatic Sciences, 74(5), 636–649. doi: 10.1139/cjfas-2016-0008

#source("my_dist.R")
cov.Exp <- function(x1, x2 = NULL, ls) {

	dist <- sqrt(my_squared_dist(x1, x2, ls))
    K <- exp(-dist)
	K
}
