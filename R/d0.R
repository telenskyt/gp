
# log likelihood given f (returns single number)
d0 <- function (gp, f, data = gp$data, hyperpar)
{
	if (gp$negLogLik.reindex2main && gp$GP_factor != "1") { 
		stopifnot(gpDataHasMainTable(data))
		# reindex from the GP_factor to main table 
		fact_idx <- paste0(gp$GP_factor, "_idx")
		f <- f[data[[1]][[fact_idx]]]		
	}
	par <- c(hyperpar[[".lik"]], list(f = f))
	res <- gp$negLogLik(data, par)
	if (!is.numeric(res) || !all(is.finite(res)))
		stop("result of user defined likelihood (negLogLik parameter to gp()) is not numeric: ", res)
	if (length(res) != 1)
		stop("length of the result of user defined likelihood (negLogLik parameter to gp()) is ", length(res), ", should be 1")
	-res
}