# identity covariance matrix (a.k.a. independent and identically distributed cov. matrix :)) for a given grouping factor(s).
# Could be also called a random effect (with respect to the given grouping factor)
#
# Since this function, unlike other cov.* functions, already considers the grouping factor(s), no reindexing will be needed.
#
# x1, x2 - grouping factors

cov.I.factor <- function (x1, x2 = NULL)
{
	if (is.null(x2))
		x2 <- x1
	outer(x1, x2, "==")*1L
}


# cov.I.d1 doesn't exist - there are no parameters to differentiate

