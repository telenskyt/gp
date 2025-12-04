# identity covariance matrix (a.k.a. independent and identically distributed cov. matrix :))
# Basically an independent noise, can be interpreted as overdispersion, or perhaps just residuals?
#
# nrows, ncols - dimensions; if ncols is not given, it should be symmetric matrix; if it is given, we assume that we are predicting for different data records

cov.I <- function (nrows, ncols = NULL)
{
	if (is.null(ncols))
		diag(nrows)
	else
		matrix(0, nrow = nrows, ncol = ncols)
}


# cov.I.d1 doesn't exist - there are no parameters to differentiate

