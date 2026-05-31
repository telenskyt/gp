# third derivative of log likelihood given f w.r.t f
# if the Hessian matrix is diagonal, it returns just a vector
# if it is non-diagonal, it returns a sparse matrix, which represents the 3-dimensional sparse matrix,
# but the first two dimensions are flattened into one

d3 <- function (gp, f, data = gp$data, hyperpar, fisher = NULL)
{
	mstart(id = "d3")
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
	
	F3 <- F$jacfun()$jacfun(sparse = TRUE)$jacfun(sparse = TRUE) # analogy of d2()
	res <- F3(f)
	if (!is_numeric_finite(res))
		stop("third derivative of the user defined likelihood (lik parameter to gp()) is not purely numeric and finite matrix")
	if (gp$W.type == "diag") {
		#mstart(id = "check isDiagonal(hessian)")
		if (!Matrix::isDiagonal(res))
			stop("third-order hessian matrix W is not diagonal - you need to set gp(W.type = ) to some other value than \"diag\"")
		#mstop(id = "check isDiagonal(hessian)", report_id = TRUE)
		res <- diag(res)
		if (length(res) != gp$GP_size)
			stop("length of the diagonal of the third derivative of the user defined likelihood (lik parameter to gp()) is ", length(res), ", should be ", gp$GP_size)
	} else {
		if (ncol(res) != gp$GP_size)
			stop("ncols of the third derivative of the user defined likelihood (lik parameter to gp()) are ", ncol(res), ", should be = ", gp$GP_size)
	}
	mstop(id = "d3", report_id = TRUE)
	-res		
}
