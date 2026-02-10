# like cov.I.factor, but with different sigma2 per factor level
#
# Since this function, unlike other cov.* functions, already considers the grouping factor(s), no reindexing will be needed.
#
# x1, x2 - grouping factors

cov.I.factor.sigma2 <- function (x1, x2 = NULL, sigma2_cat)
{
	if (is.null(x2))
		x2 <- x1
	sigma2_x1 <- sigma2_cat[as.integer(x1)]
	outer(x1, x2, "==")*1L*matrix(sigma2_x1, nrow = length(x1), ncol = length(x2))
}

cov.I.factor.sigma2.d1 <- function (x1, sigma2_cat, der_wrt, der_wrt_i)
{
	stopifnot(der_wrt == "sigma2_cat")
	stopifnot(der_wrt_i >= 1 && der_wrt_i <= length(sigma2_cat))
	x1_indicator <- as.integer(x1) == der_wrt_i
	outer(x1_indicator, x1_indicator, "&")*1L
}


# cov.I.d1 doesn't exist - there are no parameters to differentiate

