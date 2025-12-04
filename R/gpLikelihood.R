
# zatim nechavam data jako separatni parametr, i kdyz by mozna stacilo brat proste data z gp$data
# (ale co cross-validace treba?)
# hyperparametry  pro likelihood zatim vubec neresim




d0 <- function (gp, f, data = gp$data, hyperpar = numeric(0))
{
	if (gp$GP_factor != "1") { 
		# reindex from the GP_factor to main table 
		fact_idx <- paste0(gp$GP_factor, "_idx")
		f <- f[data[[1]][[fact_idx]]]		
	}
	gp$ll(data, f)
}

d1 <- function (gp, f, data = gp$data, hyperpar = numeric(0))
{
	if (gp$GP_factor != "1") { 
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
	
	drop(F$jacobian(f))
}


d2 <- function (gp, f, data = gp$data, hyperpar = numeric(0))
{
	if (gp$GP_factor != "1") { 
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
	diag(F2(f))	
}

d3 <- function (gp, f, data = gp$data, hyperpar = numeric(0))
{
	if (gp$GP_factor != "1") { 
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
	diag(F3(f))
	
}

# hyperparametry  pro likelihood zatim vubec neresim
hyperpar_first_likelihood_ind <- function(gp)
{
	return (Inf)
}