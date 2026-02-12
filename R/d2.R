# second derivative of log likelihood given f w.r.t f: at the moment, just a diagonal of Hessian matrix (marked as W in R&W 2006)
d2 <- function (gp, f, data = gp$data, hyperpar)
{
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
	
	F2 <- F$jacfun()$jacfun(sparse = TRUE) # Kasper's advice on how to calculate just Hessian diagonal! https://groups.google.com/g/tmb-users/c/fAaEhwW1niU
	hes <- F2(f)
	#mstart(id = "check isDiagonal(hessian)")
	stopifnot(Matrix::isDiagonal(hes))
	#mstop(id = "check isDiagonal(hessian)", report_id = TRUE)
	res <- diag(hes)

	if (!is.numeric(res) || !all(is.finite(res)))
		stop("second derivative of the user defined likelihood (negLogLik parameter to gp()) is not numeric: ", res)
	if (length(res) != gp$GP_size)
		stop("length of the second derivative the of user defined likelihood (negLogLik parameter to gp()) is ", length(res), ", should be ", gp$GP_size)
	-res	
}
