
#library(methods)

#' @import methods
#' xxx@importFrom methods setRefClass
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
		likFun = "ANY"
	),
	methods = list(
	
initialize = function(
	terms = function (data, par) NULL, 
	link = function (data, par) par, 
	process, 
	lik)
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
},


### template interface

# par - a list of vectors or matrices
# assumption - each vector element (or matrix row) behaves independently
lik = function(par, log = FALSE)
{
	if (!is.null(likFun))
		return(likFun(par))

	if (!is.null(likType))
		return(do.call(paste0("lik", capFirst(likType)), list(par = par, log = log)))

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
nll = function (data, par) 
{
	stages(data, par, stages = 1:5)
},

# pick which stages you want to run
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
	if (4 %in% stages)
		par <- lik(par)
	if (5 %in% stages)
		par <- -sum(log(par))
	return(par)
}


))

if (0) {
a <- likTempPhased$new(process = function () {}, lik = function () {})

a <- likTempPhased$new(process = function () {}, lik = "bern")
}
