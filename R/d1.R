
d1 <- function (gp, f, data = gp$data, hyperpar = numeric(0))
{
	if (gp$negLogLik.reindex2main && gp$GP_factor != "1") { 
		stopifnot(gpDataHasMainTable(data))
		# reindex from the GP_factor to main table in the closure wrapper, so that also the automatic differentiation
		# gives vector of the dimension of GP_factor!
		cmb <- function(func, data) function(f) {
			fact_idx <- paste0(gp$GP_factor, "_idx")
			f2 <- f[data[[1]][[fact_idx]]]
			func(data, f2)
		}
	}
	else 
		cmb <- function(func, data) function(f) func(data, f)
			
	# cmb - the closure trick somewhere from the RTMB docs, to tie this to particular data and prevent unnecessary headaches!
	#mstart(id = "MakeTape", mem_precise = TRUE)
	F <- MakeTape(cmb(gp$negLogLik, data), numeric(gp$GP_size))
	#mstop(id = "MakeTape")
	
	res <- drop(F$jacobian(f))
	if (!is.numeric(res) || !all(is.finite(res)))
		stop("first derivative of the user defined likelihood (negLogLik parameter to gp()) is not numeric: ", res)
	if (length(res) != gp$GP_size)
		stop("length of the first derivative the of user defined likelihood (negLogLik parameter to gp()) is ", length(res), ", should be ", gp$GP_size)
	-res	
}

