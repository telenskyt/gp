
# zatim nechavam data jako separatni parametr, i kdyz by mozna stacilo brat proste data z gp$data
# (ale co cross-validace treba?)
# hyperparametry  pro likelihood zatim vubec neresim




d0 <- function (gp, f, data = gp$data, hyperpar = numeric(0))
{
	if (gp$logLik.reindex2main && gp$GP_factor != "1") { 
		stopifnot(gpDataHasMainTable(data))
		# reindex from the GP_factor to main table 
		fact_idx <- paste0(gp$GP_factor, "_idx")
		f <- f[data[[1]][[fact_idx]]]		
	}
	res <- gp$ll(data, f)
	if (!is.numeric(res) || !all(is.finite(res)))
		stop("result of user defined likelihood (logLik parameter to gp()) is not numeric: ", res)
	if (length(res) != 1)
		stop("length of the result of user defined likelihood (logLik parameter to gp()) is ", length(res), ", should be 1")
	res
}

d1 <- function (gp, f, data = gp$data, hyperpar = numeric(0))
{
	if (gp$logLik.reindex2main && gp$GP_factor != "1") { 
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
	F <- MakeTape(cmb(gp$ll, data), numeric(gp$GP_size))
	#mstop(id = "MakeTape")
	
	res <- drop(F$jacobian(f))
	if (!is.numeric(res) || !all(is.finite(res)))
		stop("first derivative of the user defined likelihood (logLik parameter to gp()) is not numeric: ", res)
	if (length(res) != gp$GP_size)
		stop("length of the first derivative the of user defined likelihood (logLik parameter to gp()) is ", length(res), ", should be ", gp$GP_size)
	res	
}


d2 <- function (gp, f, data = gp$data, hyperpar = numeric(0))
{
	if (gp$logLik.reindex2main && gp$GP_factor != "1") { 
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
	F <- MakeTape(cmb(gp$ll, data), numeric(gp$GP_size))
	#mstop(id = "MakeTape")
	
	F2 <- F$jacfun()$jacfun(sparse = TRUE) # Kasper's advice on how to calculate just Hessian diagonal! https://groups.google.com/g/tmb-users/c/fAaEhwW1niU
	res <- diag(F2(f))	
	if (!is.numeric(res) || !all(is.finite(res)))
		stop("second derivative of the user defined likelihood (logLik parameter to gp()) is not numeric: ", res)
	if (length(res) != gp$GP_size)
		stop("length of the second derivative the of user defined likelihood (logLik parameter to gp()) is ", length(res), ", should be ", gp$GP_size)
	res	

}

d3 <- function (gp, f, data = gp$data, hyperpar = numeric(0))
{
	if (gp$logLik.reindex2main && gp$GP_factor != "1") { 
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
	F <- MakeTape(cmb(gp$ll, data), numeric(gp$GP_size))
	#mstop(id = "MakeTape")
	
	F3 <- F$jacfun()$jacfun(sparse = TRUE)$jacfun(sparse = TRUE) # analogy of d2()
	res <- diag(F3(f))
	if (!is.numeric(res) || !all(is.finite(res)))
		stop("third derivative of the user defined likelihood (logLik parameter to gp()) is not numeric: ", res)
	if (length(res) != gp$GP_size)
			stop("length of the third derivative the of user defined likelihood (logLik parameter to gp()) is ", length(res), ", should be ", gp$GP_size)
	res		
}

# hyperparametry  pro likelihood zatim vubec neresim
hyperpar_first_likelihood_ind <- function(gp)
{
	return (Inf)
}