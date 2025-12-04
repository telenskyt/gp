

#source("cov_funcs.R")

# ????
# K_matrix() => K()??
# dK_dli() => dK()?? -> dK_dhyp


# K_cache() - build a cache to avoid duplicate computations in multiple calls to K_matrix()/dK_dli() in specific scenarios
#
# The caller specifies the scenario using the `cache_for` argument. It is the responsibility of the caller to make sure the cache is always valid and up to date.
# 
# hyperpar - list of lists of vectors
#
# cache_for: for what use you want to use the cache?
#	"derivative" - dK_dhi() function will be called along with K_matrix(), for the exact same hyperparameters and data (??? bude tam vystupovat x2 nebo jen x1 - symetricke matice?)
#	"gradient" - K_matrix() will be called many times along a gradient of one or more variables specified by the `gradient_along` parameter, for the exact same hyperparameters and data 
#
# gradient_along - named list (names correspond to table names from the x1/x2 dataset(s)) of vectors of columns from a given table (either column names or indices).
#	This is used for cache_for = "gradient" to specify that all the other data columns except these will be fixed for the given caching session.
#
# components - see validate_components()

# !!! slo by tam mit md5 tech dat a hyperparametru pro jistotu

K_cache <- function (gp, hyperpar, x1, x2 = NULL, cache_for = c("derivative", "gradient"), gradient_along = NULL, components = NULL)
{

	# validate the parameters
	cache_for <- match.arg(cache_for)	
	stopifnot(class(x1) == "gpData")
	stopifnot(is.null(x2) || class(x2) == "gpData")
	ccs <- validate_components(gp, components) # cc - covariance component	
	
	if (cache_for == "derivative") {
		# cache all cov.SE and cov.NN.additive components
		
		K.cache <- list(
			cache_for = cache_for,
			components = list(),
			hyperpar = hyperpar # just for check; ideally, we should also check if the data are the same, but md5sum takes too much time
		)
		
		ccs <- (gpComponentsTable(gp, x1, components) %>% 
			filter(component %in% ccs, cov_fun %in% c("cov.SE", "cov.NN.additive")) %>% 
			select(component))[,1] 
					
		for (cc in ccs) { 
			cov_fun_name <- gp$covComp[[cc]]$cov_fun
			mat <- gp$covComp[[cc]]$mat
			fact <- gp$covComp[[cc]]$fact
			cov_fn_args <- hyperpar[[cc]]
			cov_fn_args[['sigma2']] <- NULL 
			x1m <- x1[[mat]]
			x2m <- if (is.null(x2)) NULL else x2[[mat]]
			K.cache$components[[cc]] <- do.call(cov_fun_name, c(list(x1 = x1m, x2 = x2m), cov_fn_args))
		}		
		return(K.cache)
	} else if (cache_for == "gradient") {
		return(NULL) # prozatimni reseni - taky moznost! Neni zatim implementovano 
	}
}


# hyperpar - list of lists of vectors
# components - see validate_components()

# DIMENZE VYSTUPU:
#
# varianta 1: stanovit dimenzi K dle dimenze datove sady, s prihlednutim k defaultni dimenzi GP (tj. zda se bude re-indexovat na uroven main table se rozhoduje dle datove sady)
# 
# - pokud v danych datech vsechny tables maji stejny grouping factor, zadna table neni main, vsechny pouzite tables have same number of rows a to bude vystupni dimenze K_matrix
# - pokud v danych datech jsou ruzne grouping factors:
# 	- pokud dimenze main GP (s training dataset a vsemi komponentami) je per factor, tak vystupni dimenze bude odpovidat tomuto faktoru v dane datove sade
# 	- v opacnem pripade vystupni dimenze bude odpovidat main table dane datove sady
# 
# => tj. bude to jednoduse gpDataSize(x, gp$GP_factor) !! Kde x je dana datova sada
# 
# 
# ! bacha na specialni pripad, kdy budu delat K_matrix jen pro intercept - ale to bude asi stejne?
# 
# 
# varianta 2: stanovit dimenzi K_matrix nikoli dle dimenze datove sady, ale jen podle specifikovanych komponent (tj. ne-reindexuje se, kdyz vsechny specifikovane komponenty patri 
# jednomu grouping factoru)
# 
# - toto mi bylo sympaticke, protoze by to stanoveni dimenze asi bylo jednodussi, ale bojim se, aby to neco nerozbilo (tak radeji varianta 1, ktera zachovava puvodni funkcnost)
# 	- ono totiz kdyz bych si pak delal K_matrix pro jednotlive komponenty, tak by to melo ruzne dimenze, a to nevim zda by neco nerozbilo

K_matrix <- function (gp, hyperpar = gpHyperparList(gp), x1 = gp$data, x2 = NULL, K.cache = NULL, components = NULL, comp_means = FALSE, comp_means_weights = NULL)
{
	# validate the parameters
	stopifnot(class(x1) == "gpData")
	stopifnot(is.null(x2) || class(x2) == "gpData")
	ccs <- validate_components(gp, components) # cc - covariance component	
	if (!is.null(K.cache)) {
		stopifnot(K.cache$cache_for == "derivative") # "gradient" not implemented yet
		stopifnot(identical(hyperpar, K.cache$hyperpar))
	}
	# extra checks to let the user know (maybe should be just warnings?)
	if (is.null(attr(x1, "gpDataPrepared")))
		stop("perhaps you want to run gpDataPrepare() first for the x1 argument?")
	if (!is.null(x2) && is.null(attr(x2, "gpDataPrepared")))
		stop("perhaps you want to run gpDataPrepare() first for the x1 argument?")
	
	# arrange the components from smallest to largest, and the intercept (which has nrow = NA) will go last!
	ccs <- (gpComponentsTable(gp, x1, components) %>% 
		filter(component %in% ccs) %>% 
		arrange(nrow) %>% # we use nrow just from x1, not x2, but it's not that important
			# !!!! tady to razeni melo byt kvuli optimalizaci, ale ted je to jedno, protoze tam ta optimalizace neni. Ale mohla by se udelat!
			# ze by se nereindexovalo kazde Kc zvlast,ale nejdriv by se poscitaly ty mensi matice bez reindexace a reindexovalo by se to na GP scale
			# az kdyby dalsi matice byla nejaka vetsi (v tom pripade to razeni by nemelo bejt podle poctu radek, ale tak, aby se reindexovalo co nejpozdeji)
		select(component))[,1] 
		
	if (comp_means && !is.null(x2))
		stop("comp_means has only meaning for the training covariance matrix")
	
	# output matrix dimensions
	nrows <- gpDataSize(x1, fact = gp$GP_factor)
	ncols <- if (is.null(x2)) NULL else gpDataSize(x2, fact = gp$GP_factor)
		
	K <- NULL
	meansm <- NULL
	for (cc in ccs) { 
		#cat("Doing component", cc, " \n")
		cov_fun_name <- gp$covComp[[cc]]$cov_fun
		fact <- gp$covComp[[cc]]$fact

		reindex_needed <- TRUE
		if (!is.na(fact) && fact == gp$GP_factor) # we are at the scale of GP, no reindex needed
			reindex_needed <- FALSE
		
		# check - toto by melo ted platit, diky nasim omezenim:
		if (!is.na(fact) && fact != "1" && gp$GP_factor != "1")
			stopifnot(fact == gp$GP_factor)
		
		if (!is.null(K.cache) && cc %in% names(K.cache$components)) {
			Kc <- hyperpar[[cc]]$sigma2 * K.cache$components[[cc]]
		} else if (cov_fun_name == "1") { # intercept
			Kc <- matrix(hyperpar[[cc]]$sigma2, nrow = nrows, ncol = if (is.null(ncols)) nrows else ncols)
				# mrzi me, ze je tohle potreba delat pro hloupej intercept, ale mel jsem to tak i v puvodni knihovne (v K-tsc_spat2.R), nejspis asi jen 
				# kvuli tem calc_K_mean()			
			#cat("intercept dim:\n")
			#print(dim(Kc))
			reindex_needed <- FALSE
		#} else if (cov_fun_name == "I") {
		#	stopifnot(!is.na(fact) && !is.null(fact)) # should only be 1 or some true factor
		#	if (fact == "1" 
		} else if (cov_fun_name == "cov.I") { # I
			Kc <- hyperpar[[cc]]$sigma2 * cov.I(nrows, ncols) # I|fact
			reindex_needed <- FALSE
		} else if (cov_fun_name == "cov.I.factor") {
			stopifnot(gp$GP_factor == "1") # jiny pripad neni zatim osetreny, jsou tam challenges, viz GP_factor-k_rozreseni.docx
			Kc <- hyperpar[[cc]]$sigma2 * cov.I.factor(x1[[1]][[fact]], if (is.null(x2)) NULL else x2[[1]][[fact]])
			reindex_needed <- FALSE
		} else {
			mat <- gp$covComp[[cc]]$mat
			#cov_funcs[[cov_fun_name]]
			#cov_fun <- get(cov_fun_name)
			cov_fn_args <- hyperpar[[cc]]
			cov_fn_args[['sigma2']] <- NULL 
			stopifnot(!is.na(mat))
			x1m <- x1[[mat]]
			x2m <- if (is.null(x2)) NULL else x2[[mat]]
			Kc <- hyperpar[[cc]]$sigma2 * do.call(cov_fun_name, c(list(x1 = x1m, x2 = x2m), cov_fn_args))
		}
		# v nasledujicim (ale mozna i v tom vyse) bude potreba specialne osetrit intercept a cov.iid() - proste ty, co maj `mat` = NA! Nebo ne?
		if (comp_means)
			meansm_c <- calc_K_mean(Kc, cc, comp_means_weights)
		if (!reindex_needed) { # no reindex needed
			if (is.null(K))
				K <- Kc
			else
				K <- K + Kc
		} else { # we will need to re-index from `fact` to the main table
			fact_idx <- paste0(fact, "_idx")
			if (is.null(x2)) {
				if (is.null(K))
					K <- Kc[x1[[1]][[fact_idx]], x1[[1]][[fact_idx]]]
				else
					K <- K + Kc[x1[[1]][[fact_idx]], x1[[1]][[fact_idx]]]
			} else {
				if (is.null(K))
					K <- Kc[x1[[1]][[fact_idx]], x2[[1]][[fact_idx]]]
				else
					K <- K + Kc[x1[[1]][[fact_idx]], x2[[1]][[fact_idx]]]
			}
			if (comp_means)
				meansm_c <- meansm_c[x1[[1]][[fact_idx]],] # re-index rows from `fact` to the main table
		}
		#cat("component dimension: ")
		#print(dim(Kc))
		if (comp_means)
			meansm <- cbind(meansm, meansm_c)
			
		rm(Kc) # 2025-10: tento pristup bude zrat min pameti, ale melo by se nekde otestovat, jestli to pravidelne gc() neni casove moc narocne
		gc()
	}
	if (comp_means) {	
		attr(K, 'comp_means') <- meansm
	}
	K	
}


# internal, auxilliary function
calc_K_mean <- function (k, component, weights = NULL) 
{
	if (!is.matrix(k))
		return(NULL)
	else if (is.null(weights)) {
		#m <- matrix(apply(k, 1, mean), ncol = 1)
		m <- matrix(rowMeans(k), ncol = 1) # comp_means calculation took 0.3/0.2 sec 
	} else if (1) {
		stopifnot(length(weights) == nrow(k))
		m <- k %*% (weights / sum(weights))  # comp_means calculation took 0.2 sec
	} else {
		stopifnot(length(weights) == nrow(k))
		m <- matrix(rowMeans(k * matrix(weights/sum(weights), nrow = nrow(k), ncol = ncol(k), byrow = TRUE)), ncol = 1)
			#  comp_means calculation took 1.8 sec
	}
	colnames(m) <- component
	#cat("calc_K_Mean:\n")
	#print(head(m))
	m
}


# dK/dhyperparameter[i]: derivation of covariance matrix w.r.t hyperparameter with index i (within the optimized ones)
dK_dhi <- function (gp, hyperpar, x1, i, K.cache = NULL)
{
	stopifnot(class(x1) == "gpData")
	
	i <- gpHyperparIdx(gp, i)
	
	der_wrt <- gp$hyperpar[i,] # which hyperpar we differentiate w.r.t to :-)
	
	cc <- der_wrt$component
	cov_fun_name <- gp$covComp[[cc]]$cov_fun # nikoli der_wrt$cov_fun, to je NA pro sigma2!
	cov_fun_name_d1 <- paste0(der_wrt$cov_fun, ".d1")
	fact <- gp$covComp[[cc]]$fact
	
	mat <- gp$covComp[[cc]]$mat
	cov_fn_args <- hyperpar[[cc]]
	cov_fn_args[['sigma2']] <- NULL 	
	
	N <- gpDataSize(x1, fact = gp$GP_factor)
	
	reindex_needed <- TRUE
	if (!is.na(fact) && fact == gp$GP_factor) # we are at the scale of GP, no reindex needed
		reindex_needed <- FALSE
		
	if (der_wrt$hyperpar == "sigma2") {
		# just give the normal covariance matrix, but not multiplied by this sigma2
		if (!is.null(K.cache) && cc %in% names(K.cache$components)) {
			stopifnot(identical(hyperpar, K.cache$hyperpar))
			K <- K.cache$components[[cc]]
		} else if (cov_fun_name == "1") {
			K <- matrix(1, nrow = N, ncol = N)
			reindex_needed <- FALSE
		} else if (cov_fun_name == "cov.I") {
			K <- cov.I(N)
			reindex_needed <- FALSE
		} else if (cov_fun_name == "cov.I.factor") {
			stopifnot(gp$GP_factor == "1") # jiny pripad neni zatim osetreny, jsou tam challenges, viz GP_factor-k_rozreseni.docx		
			K <- cov.I.factor(x1[[1]][[fact]])
			reindex_needed <- FALSE
		} else {
			stopifnot(!is.na(mat))
			x1m <- x1[[mat]]
			K <- do.call(cov_fun_name, c(list(x1 = x1m, x2 = NULL), cov_fn_args))
		}
	} else if (cov_fun_name == "cov.SE") { # this one needs special treatment
		stopifnot(!is.null(K.cache) && cc %in% names(K.cache$components))
		stopifnot(identical(hyperpar, K.cache$hyperpar))
		stopifnot(!is.na(mat))
		x1m <- x1[[mat]]		
		K <- hyperpar[[cc]]$sigma2 * cov.SE.d1(K = K.cache$components[[cc]], x1 = x1m, ls = hyperpar[[cc]]$ls, der_wrt_i = der_wrt$i)
	} else { # derivative of generic cov.* functions w.r.t. one their internal hyperparameter
		stopifnot(!cov_fun_name %in% c("1", "cov.I", "cov.I.factor"))
		stopifnot(!is.na(mat))
		x1m <- x1[[mat]]
		K <- hyperpar[[cc]]$sigma2 * do.call(cov_fun_name_d1, c(list(x1 = x1m, der_wrt = der_wrt$hyperpar, der_wrt_i = der_wrt$i), cov_fn_args))
	}
	if (reindex_needed) { # we will need to re-index from `fact` to the main table
		fact_idx <- paste0(fact, "_idx")
		K <- K[x1[[1]][[fact_idx]], x1[[1]][[fact_idx]]]
	}
	K
}


