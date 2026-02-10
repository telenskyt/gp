
#source("cov_funcs.R")

# f - formula

gpFormula <- function (f)
{
	# proceeding according to https://stackoverflow.com/a/65720737/

	if (f[[1]] != "~")
		stop("formula must contain ~")

	if (length(f) != 2) 
		stop("left side of the formula should be empty")


	f <- f[[2]]

	covComp <- list()  # covariance components
	dataReq <- list() # data requirements

	repeat { # over summed terms

	if (f[[1]] == '+') {
		stopifnot(length(f) == 3)
		t <- f[[3]]
		f <- f[[2]]
	} else {
		t <- f
		f <- NULL
	}
	# t is the term now

	if (t[[1]] != ':') {
		stop("Expecting ':' in the term '", deparse(t), "'")
	}
	if (length(t) != 3) {
		stop("The term '", deparse(t), "' should have right and left side around ':'")
	}

	name <- deparse(t[[2]]) # component name
	if (name %in% names(covComp)) {
		stop("Component name '", name, "' used twice. Can't use the same name for multiple components.")
	}

	if (deparse(t[[3]]) == "1") { # intercept - <name>:1
		cc <- list()
		cc[[name]] <- list(
			cov_fun = "1",
			mat = NA,
			fact = NA # here it is a special case - it is not factor "1", which means matrix without grouping factor, because it is basically dimensionless!
		)
	} else if (deparse(t[[3]]) == "I") { # <name>:I
		cc <- list()
		cc[[name]] <- list(
			cov_fun = "cov.I",
			mat = NA,
			fact = "1"
		)
	} else if (length(t[[3]]) == 1) { # <name>:fact
		fact <- deparse(t[[3]])
		cc <- list()
		cc[[name]] <- list(
			cov_fun = "cov.I.factor.sigma2",
			mat = NA,
			fact = fact
		)		
	} else if (t[[c(3,1)]] == '(') { # <name>:(I|fact)
		stopifnot(length(t[[c(3,2)]]) == 3)
		if (t[[c(3,2,1)]] != '|') 
			stop("Symbol '|' expected in the term '", deparse(t[[c(3,2)]]), "'")
		if (deparse(t[[c(3,2,2)]]) != "I") 
			stop("Expected 'I' before '|' in the term '", deparse(t[[c(3,2)]]), "'")
		fact <- deparse(t[[c(3,2,3)]])
		cc <- list()
		cc[[name]] <- list(
			cov_fun = "cov.I.factor",
			mat = NA,
			fact = fact
		)
	} else if (length(t[[3]]) == 2 && deparse(t[[c(3,1)]]) %in% names(cov_funcs)) { # <name>:cov.fn(par)
		cov_fun <- deparse(t[[c(3,1)]])
		scaling <- cov_funcs[[cov_fun]]$scaling
		if (length(t[[c(3,2)]]) == 1) { # `par` is just matrix name without |
			mat <- deparse(t[[c(3,2)]])
			fact <- "1"
		} else { # `par` is mat|fact
			if (length(t[[c(3,2)]]) != 3)
				stop("Parse error in the term '", deparse(t[[c(3,2)]]), "'")
			if (t[[c(3,2,1)]] != '|') 
				stop("Symbol '|' expected in the term '", deparse(t[[c(3,2)]]), "'")
			mat <- deparse(t[[c(3,2,2)]])
			fact <- deparse(t[[c(3,2,3)]])
		}
		cc <- list()
		cc[[name]] <- list(
			cov_fun = cov_fun,
			mat = mat,
			fact = fact
		)
		if (!mat %in% names(dataReq$mats)) { # vloz matrix s faktorem fact a zkontroluj ze tam neni s zadnym jinym
			dr <- list()
			dr[[mat]] <- list(fact = fact, scaling = scaling)
			dataReq$mats <- c(dr, dataReq$mats) # add before the existing list, so that order is preserved
		} else {
			stopifnot(!is.na(dataReq$mats[[mat]]))
			if (!compareNA(dataReq$mats[[mat]]$fact, fact))
				stop("Table '", mat, "' cannot be used with two different factors")
			if (dataReq$mats[[mat]]$scaling != scaling)
				stop("Table '", mat, "' is used in two different covariance functions, one requires scaling = TRUE and one scaling = FALSE. This is not supported at the moment - you need to provide the data in two separate tables (duplicate it).")
		}
	} else {
		stop("Failed to parse the term '", deparse(t[[3]]), "'")
	}
	
	covComp <- c(cc, covComp) # add before the existing list, so that order is preserved	
	fact <- cc[[name]]$fact
	if (!is.na(fact) && fact != "1")
		dataReq$factors <- union(dataReq$factors, fact)
		
	if (is.null(f))
		break

	}

	x <- list(
		covComp = covComp,
		dataReq = dataReq
	)
	class(x) <- "gp"
	return(x)
}



if (0) {

f0
sxp(covComp)
lobstr::tree(covComp)

str(covComp)
str(dataReq)

}
