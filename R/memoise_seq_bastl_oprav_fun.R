# sequential version of memoise()
# bastl: pribastlil jsem si tam specialitu ciste pro muj soucasnej projekt :-D
# (hledej string bastl)
# druha verze: oprav_fun - jsou-li mezi argumenty funkce, obcas to stejny call nerozpozna jako stejny protoze se lisi jejich environmenty...
#              tak tam misto cele fce i s env dam do serializace jen jejich args a body (snaha o prime odstraneni env. nezabrala)
#
#
# 2025-11 changes:
#	- switched from base::body() to rlang::fn_body(), so that it ignores the srcref attributes when options(keep.source = TRUE)
#	- switched from comparing serialize() outputs to use identical(), because the ignore.bytecode = TRUE is essential! There were differences in bytecode for some reason
#	  and so the memoise'd result wasn't used.

#library(pryr) # make_call; driv nebylo potreba toto sem davat
#library(rlang) # fn_body()
#	!!! it is important to use rlang::fn_body() instead of body(), because body() would include all the srcref attributes when options(keep.source = TRUE)!!!

#'@importFrom pryr make_call
#'@importFrom rlang fn_body
memoise_seq_bastl <- function (f, ..., envir = environment(f), cache = cache_memory()) 
{
	
	cat("memoise ver 6\n")
	f_formals <- formals(args(f))
	if (is.memoised(f)) {
		stop("`f` must not be memoised.", call. = FALSE)
	}
	f_formal_names <- names(f_formals)
	f_formal_name_list <- lapply(f_formal_names, as.name)
	init_call_args <- setNames(f_formal_name_list, f_formal_names)
	init_call <- pryr::make_call(quote(`_f`), init_call_args)
	validate_formulas(...)
	additional <- list(...)
	memo_f <- eval(bquote(function(...) {
		called_args <- as.list(match.call())[-1]
		default_args <- Filter(function(x) !identical(x, quote(expr = )), 
							   as.list(formals()))
		default_args <- default_args[setdiff(names(default_args), 
											 names(called_args))]
		args <- c(lapply(called_args, eval, parent.frame()), 
				  lapply(default_args, eval, envir = environment()))
		call_serialized <- list(
			body = rlang::fn_body(`_f`) 
			,args = args
			,lapply(`_additional`, function(x) eval(x[[2L]], environment(x)))
		)
																		   
		#cat("in memoise'd fun: counter =", `_counter`, "\n")
		#`_counter` <- `_counter` + 1
		#pruchod <<- pruchod + 1
		#if (pruchod > 1) {
		#ahoj <- nazdar
		#}
		# trochu sileny zpusob pouziti cache, ale nevim, jak jinak na nej!
		# ono totiz nejde menit ty promenne v environmentu!! (Jak jsem delal ten pokus s tim `_counter`, zustaval stale na nule!)
		# Aha! uz asi vim, proc to neslo! Ja blbec zapomnel pouzit <<- misto <-! To by asi fungovalo! ALe tak necham uz to takhle,
		# pres tu cache. Ono to ten environment udrzi malej aspon.
		if (`_cache`$has_key("_MYDATA_"))
			data <- `_cache`$get("_MYDATA_")
		else
			data <- list(`_last_call` = "")
		if (0) { # debug
			body <- rlang::fn_body(`_f`)
			addit <- lapply(`_additional`, function(x) eval(x[[2L]], 
																			   environment(x)))
			save(body, args, addit, file = paste0("tmp_iter", hyper_iter, ".Rdata"))		
		
			sink(file = paste0("tmp_iter", hyper_iter, "-ser.txt"))
			str(unserialize(call_serialized))
			sink()
			sink(file = paste0("tmp_iter", hyper_iter, "-last.txt"))
			if (length(data$`_last_call`) > 1)
				str(unserialize(data$`_last_call`))
			sink()
		}
		if (identical(call_serialized, data$`_last_call`)) { # porovnavam cely objekt, nejen hash. Tak mam 100% jistotu.
			res <- data$`_last_res`
			cat("memoise: using memoised res. Environment size =", object.size(environment())/1024, "kB\n")
		}
		else {
			#cat("** BEFORE memoise reset:\n")
			#print(gc())
			data$`_last_res` <- NULL # dulezite aby se uvolnila pamet!
			`_cache`$set("_MYDATA_", NULL) # dulezite! Nejdriv smazat, protoze data jsou v mem pripade velka! Slo by i volat $reset()
											# gc() at si zavola ta funkce sama kdyz bude potrebovat
			#`_cache`$reset()
			#cat("\nresetting the shit\n")
			#cat("** AFTER memoise reset:\n")
			#print(gc())
		#if (hyper_iter > 0) {
		#ahoj <- nazdar
		#}			
			cat("memoise: calling the function\n")
			#browser()
			#if (hyper_iter > 0)
			#	recover()
			res <- withVisible(.(init_call))
			data$`_last_res` <- res
			data$`_last_call` <- call_serialized
			`_cache`$set("_MYDATA_", data)
			{ # bastlim!!!
				hyper_iter <<- hyper_iter + 1
			}
		}
		if (res$visible) {
			res$value
		}
		else {
			invisible(res$value)
		}
	}, as.environment(list(init_call = init_call))))
	formals(memo_f) <- f_formals
	attr(memo_f, "memoised") <- TRUE
	if (is.null(envir)) {
		envir <- baseenv()
	}
	memo_f_env <- new.env(parent = envir)
	memo_f_env$`_cache` <- cache
	#	memo_f_env$`_last_call` <- ""
	#	memo_f_env$`_last_res` <- NA
	#	memo_f_env$`_counter` <- 0
	memo_f_env$`_f` <- f
	memo_f_env$`_additional` <- additional
	environment(memo_f) <- memo_f_env
	class(memo_f) <- c("memoised", "function")
	memo_f
}







