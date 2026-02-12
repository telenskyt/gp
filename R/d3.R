# third derivative of log likelihood given f w.r.t f: at the moment, since the Hessian matrix is diagonal (marked as W in R&W 2006), this is just a vector
d3 <- function (gp, f, data = gp$data, hyperpar)
{
	mstart(id = "d3")
	cmb <- function(func, data) function(f) {
		if (gp$negLogLik.reindex2main && gp$GP_factor != "1") { 
			stopifnot(gpDataHasMainTable(data))	
			# reindex from the GP_factor to main table in this closure wrapper, so that also the automatic differentiation
			# gives vector of the dimension of GP_factor!
			fact_idx <- paste0(gp$GP_factor, "_idx")
			f <- f[data[[1]][[fact_idx]]]
		}
		par <- c(hyperpar[[".lik"]], list(f = f))
		func(data, par)
	}			
	# cmb - the closure trick somewhere from the RTMB docs, to tie this to particular data and prevent unnecessary headaches!
	#mstart(id = "MakeTape", mem_precise = TRUE)
	F <- MakeTape(cmb(gp$negLogLik, data), numeric(gp$GP_size))
	#mstop(id = "MakeTape")
	
	F3 <- F$jacfun()$jacfun(sparse = TRUE)$jacfun(sparse = TRUE) # analogy of d2()
	res <- diag(F3(f))
	if (!is.numeric(res) || !all(is.finite(res)))
		stop("third derivative of the user defined likelihood (negLogLik parameter to gp()) is not numeric: ", res)
	if (length(res) != gp$GP_size)
			stop("length of the third derivative the of user defined likelihood (negLogLik parameter to gp()) is ", length(res), ", should be ", gp$GP_size)
	mstop(id = "d3", report_id = TRUE)
	-res		
}
