# second derivative of log likelihood given f w.r.t f: at the moment, just a diagonal of Hessian matrix (marked as W in R&W 2006)
d2 <- function (gp, f, data = gp$data, hyperpar, fisher = NULL)
{
	#mstart(id = "d2")
	if (is.null(fisher))
		nll_fun <- gp$lik$nll
	else
		nll_fun <- Fisher_dispatcher(gp, f, data, hyperpar, fisher)
	
	#mstart(id = "MakeTape", mem_precise = TRUE)
	F <- MakeTape(function(f) {
		par <- gpGetParForLikTemplate(gp, f, data, hyperpar)
		nll_fun(data, par)
	}, numeric(gp$GP_size))
	#mstop(id = "MakeTape", report_id = TRUE)

	#mstart(id = "optimizations")
	# tyto vyrazne pomuzou! (nebo aspon jedna z nich):
	if (isTRUE(attr(nll_fun, "eliminate"))) {
		#mstart(id = "eliminate")
		F$simplify("eliminate")
		#mstop(id = "eliminate", report_id = TRUE)	
	}
	if (isTRUE(attr(nll_fun, "optimize"))) {
		#mstart(id = "optimize")
		F$simplify("optimize")
		#mstop(id = "optimize", report_id = TRUE)	
	}
	#mstop(id = "optimizations", report_id = TRUE)	
	
	#mstart(id = "jacfun jacfun", mem_precise = TRUE)
	F2 <- F$jacfun()$jacfun(sparse = TRUE) # Kasper's advice on how to calculate just Hessian diagonal! https://groups.google.com/g/tmb-users/c/fAaEhwW1niU
	#mstop(id = "jacfun jacfun", report_id = TRUE)
	#mstart(id = "tape eval")
	res <- F2(f)
	#mstop(id = "tape eval", report_id = TRUE)
	if (!is_numeric_finite(res))
		stop("second derivative of the user defined likelihood (lik parameter to gp()) is not purely numeric and finite matrix")
	if (gp$W.type == "diag") {
		#mstart(id = "check isDiagonal(hessian)")
		if (!Matrix::isDiagonal(res))
			stop("hessian matrix W is not diagonal - you need to set gp(W.type = ) to some other value than \"diag\"")
		#mstop(id = "check isDiagonal(hessian)", report_id = TRUE)
		res <- diag(res)
		if (length(res) != gp$GP_size)
			stop("length of the diagonal of the second derivative of the user defined likelihood (lik parameter to gp()) is ", length(res), ", should be ", gp$GP_size)
	} else {
		if (nrow(res) != gp$GP_size || ncol(res) != gp$GP_size)
			stop("dimensions of the second derivative of the user defined likelihood (lik parameter to gp()) are ", dim(res), ", should be both = ", gp$GP_size)
	}
	#mstop(id = "d2", report_id = TRUE)
	-res	
}
