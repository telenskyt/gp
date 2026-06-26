
#library(methods)

#' @import methods
#' @export likTempPhased
#' @exportClass likTempPhased
likTempPhased <- setRefClass(
	"likTempPhased",
	#contains = "likTemp",
	fields = list(
		terms = "function",
		link = "function",
		process = "function",
		likType = "ANY",
		likFun = "ANY",
		f_start = "ANY"
	),
	methods = list(
	
initialize = function(
	terms = function (data, par) NULL, 
	link = function (data, par) par, 
	process, 
	lik,
	f_start = NULL)
{
	#if (!is.null(terms))
	terms <<- terms
	link <<- link
	process <<- process
	likType <<- NULL
	likFun <<- NULL
	if (is.character(lik))
		likType <<- lik
	else if (is.function(lik)) 
		likFun <<- lik
	if (!is.null(f_start)) {
		if (!(is.function(f_start) && identical(names(formals(f_start)), c("prev_terms", "current_terms"))))
			stop("f_start must be a function(prev_terms, current_terms)")
	}
	f_start <<- f_start
},


### template interface
#
# par - a list of vectors or matrices, contains the data as well as the parameters, in the format required for this given likelihood
# optional parameters:
#	log		return logarithm? Must be FALSE by default
#	sum		sum the logarithms? Must be FALSE by default
#
# returns: likelihood, as a probability, kept elementwise! 
#
# assumption - each vector element (or matrix row) behaves independently
lik = function(par, ...)
{
	lik.getfun()(par, ...)
},

lik.getfun = function ()
{
	if (!is.null(likFun))
		return(likFun)

	if (!is.null(likType))
		return(match.fun(paste0("lik", capFirst(likType))))
		
	stop("error")
},

# par - a list of vectors or matrices; same format as for the lik() method
# samples - number of samples
# returns the simulated data as a list of matrix(es), columns corresponding to samples
simulate = function (par, samples)
{
	if (!is.null(likType))
		return(do.call(paste0("sim", capFirst(likType)), list(par = par, samples = samples)))
	else
		stop("cannot simulate if likelihood is not specified as character string")
},

# nll: runs all the template phases one by one, and returns the negative log likelihood
#
# data: scaled data, always contains "y" - the training data
# par: 
#	- contains likelihood hyperparameters, if there are some
#	- contains "f"
#
# returns: negative log likelihood (a single number)
nll = function (data, par) 
{
	stages(data, par, stages = 1:5)
},

# pick which stages you want to run; data & par must correspond to the stages selected, as will the return value
stages = function (data, par, stages = 1:5)
{
	if (1 %in% stages) {
		# first stage, "terms", is called without "f"; the terms() must never use "f", that must be kept for further stages
		par1 <- par
		par1$f <- NULL
		par2 <- terms(data, par1) # might return NULL
		# we add "f" back and create what we can call "predictions"; .lik hyperparameters and everything else is no longer there
		par <- bind_cols(par2, data.frame(f = par$f)) # with bind_cols(), par2 can be NULL - this is important
		# par are "predictions" now 
	}
	if (2 %in% stages)
		par <- link(data, par)
	if (3 %in% stages)
		par <- process(data, par)
	if (all(4:5 %in% stages) && "log" %in% names(formals(lik.getfun()))) {
		# both stages to be run; then it is better to calculate the log likelihood directly, if possible (for numerical reasons)
		if ("sum" %in% names(formals(lik.getfun())))
			par <- -lik(par, log = TRUE, sum = TRUE) # in some cases, even summing directly might optimize things (especially RTMB::MakeTape), see e.g. likBern
		else
			par <- -sum(lik(par, log = TRUE))
	} else {
		if (4 %in% stages)
			par <- lik(par)
		if (5 %in% stages)
			par <- -sum(log(par))
	}
	return(par)
}


))

if (0) {
a <- likTempPhased$new(process = function () {}, lik = function () {})

a <- likTempPhased$new(process = function () {}, lik = "bern")
}
