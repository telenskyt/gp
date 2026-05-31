
# The functions d[0-3]* are now in separate files
#
# obecna note k temto fcim - zatim nechavam data jako separatni parametr, i kdyz by mozna stacilo brat proste data z gp$data
# (ale co cross-validace treba?)

gpGetParForLikTemplate <- function(gp, f, data, hyperpar)
{
	if (gp$lik.reindex2main && gp$GP_factor != "1") { 
		stopifnot(gpDataHasMainTable(data))	
		# reindex from the GP_factor to main table, as requested by the template
		# this should be called within the MakeTape template, so that also the automatic differentiation
		# gives vector of the dimension of GP_factor!
		fact_idx <- paste0(gp$GP_factor, "_idx")
		f <- f[data[[1]][[fact_idx]]]
	}
	par <- c(hyperpar[[".lik"]], list(f = f))
	par
}