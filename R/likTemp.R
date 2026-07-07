
#library(methods)

#' @import methods
#' @export likTempPhased
#' @exportClass likTempPhased

#' @include suppress_warnings.R
likTempPhased <- suppress_warnings_from(setRefClass(
	"likTempPhased",
	#contains = "likTemp",
	fields = list(
		terms = "function",
		link = "function",
		process = "function",
		likType = "ANY",
		likFun = "ANY",
		f_start = "ANY",
		constr.args = "list"
	),
	methods = list(
	
initialize = function(
	terms = function (data, par) NULL, 
	link = function (data, par) par, 
	process, 
	lik,
	f_start = NULL)
{
	# now, the user-supplied template functions can no longer use variables from closure! All they can use must be from data & par arguments 
	# or base functions. We are doing this to save potential space when storing the object on the disk. 
	#	- another possible solution would be to allow the closures, but warn if they are big	
	check_fun_env_size(terms)
	check_fun_env_size(link)
	check_fun_env_size(process)
	if (is.function(lik))
		check_fun_env_size(lik)
	if (is.function(f_start))
		check_fun_env_size(f_start)

	constr.args <<- as.list(environment()) # save the arguments into a list
	
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

	if (!is.null(likType)) {
		name <- paste0("lik", capFirst(likType))
		fun <- match.fun(name)
		if (!is.function(fun))
			stop("Likelihood function ", name, " for likelihood ", likType, " not found")
		return(fun)
	}
		
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
},

# Template method to evaluate predictive performance for given likelihood
#	- now an internal helper to call the predPerf function for likelihood given by likType
#
# par, par.null: `par` from the "process" stage (stage 3); par.null - dtto but from the null model
predPerfFun = function (par, par.null)
{
	if (!is.null(likType)) {
		name <- paste0("predPerf", capFirst(likType))
		fun <- match.fun(name)
		if (!is.function(fun))
			stop("Predictive performance function ", name, " for likelihood ", likType, " not found")
	} else {
		stop("predictive performance function not available when likelihood is not specified by likType") 
	}
	fun(par, par.null)
},

# This method is the interface provided for outside:
#
# Evaluate predictive performance for given predictions
#
# @param data 
# @param pred predictions to be evaluated
# @param pred.null predictions of a null model (some performance statistics need it as a baseline; see eg. ?predPerfBern)
# @details Both \code{pred} and \code{pred.null} are expected to be predictions from predict(type = "terms").
# @returns the predicted performance statistics, see eg. ?predPerfBern
# @seealso [predPerfBern()]
#
# @export
predPerf = function(data, pred, pred.null)
{
	par <- stages(data = data, par = pred, stages = 2:3)
	par.null <- stages(data = data, par = pred.null, stages = 2:3)
	predPerfFun(par, par.null)
}


)), message = "local assignment to field name will not change the field", fun = ".checkFieldsInMethod")

if (0) {
a <- likTempPhased$new(process = function () {}, lik = function () {})

a <- likTempPhased$new(process = function () {}, lik = "bern")
}


# Checks the size of the passed function along with its environment
#' @importFrom lobstr obj_size
check_fun_env_size <- function(f, max_size = getOption("gp.max_fun_size_warn", 1 * 1024^2))
{
	name <- deparse1(substitute(f))
	size <- lobstr::obj_size(f)

	if (size > max_size)
		warning(
			"Function environment of `", name, "` is large: ", format(size, units = "auto"), ". ",
			"This may make packed gp objects large. Consider defining the function ",
			"in a smaller environment. Note that all important things should be passed ",
			"through the function arguments!",
			call. = FALSE
		)

}
