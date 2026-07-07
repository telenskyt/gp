#' Perform an N-fold cross-validation of a Gaussian process model
#'
#' Function will fit the GP model for each of the N cross-validation folds, and evaluate cross-validation prediction for the whole training data set.
#'
#' @param gp GP model object
#' @param fold.col either a name of a column in the main table, or a vector along factor \code{fold.fact}.
#'		The column (or the supplied vector) must be a vector of integers from \code{1} to \code{N} (\code{N} being the number of folds), 
#'		specifying the number of the cross-validation fold the given record belongs to.
#' @param fold.fact a factor along which the cross-validation folds (\code{fold.col}) are specified. Factor \code{"1"} means the folds are specified for the rows of the main table.
#'		Note that if the GP dimension is given by some real grouping factor (i.e. \code{gp$GP_factor != "1"}), then the \code{fold.fact} must be that factor.
#' @param folds integer vector of folds to fit; if \code{NULL} (the default), all folds are fit
#' @param start.from.model an object of class \code{gp} to take the starting values from (from the \code{value} column in the hyperparameter table, 
#'		see \code{\link{gpHyperparStartFromModel}}). If specified, then for each cross-validation fold model to be fit, 
#' 		the starting values of the hyperparameters will be taken from the corresponding fold model wherever possible. This implies that the object specified in this argument
#'		must have already been fitted by gpFitCV with the exact same folds.
#' @param parallel should the cross-validation models be run in parallel? Uses foreach(), and requires parallel background for foreach() to be already registered.
#' @param fn.prefix character; prefix of file names for log files (\code{log.fn}) and dump files (\code{dump.fn}), in case these are not \code{NULL}
#' @param log.fn if not \code{NULL}, standard and error output of each fold model job will be saved into a log file with this file name. Special \code{\%} sequences can be used, see Details below. 
#' @param dump.fn if not \code{NULL}, debug dump of given fold model will be saved upon an error, with this file name (the .rda extension will be added to it). Special \code{\%} sequences can be used, see Details below. 
#' @param tr.max.lines the \code{max.lines} parameter for \code{traceback}, i.e. the maximum number of lines printed per call when error occurs
#' @param pred.options list of options to be passed to the \code{\link{predict.gp}} method.
#' @param lmFit.options list of options to be passed to \code{\link{lmFit}} when fitting the intercept-only null model for each fold.
#' @param lmFit.pred.options list of options to be passed to the \code{predict.lmFit} method when predicting from the null model.
#' @param ... options to be passed to the \code{\link{gpFit}} method.
#'
#' @details The arguments \code{log.fn} and \code{dump.fn} allow for special sequences:
#' - \code{\%f} - number of the cross-validation fold
#' - \code{\%h} - hostname, i.e. the name of the machine where the worker job runs
#' - \code{\%p} - process ID of the worker job
#'
#' @export
#' @importFrom foreach foreach getDoParRegistered getDoParName
#' @importFrom foreach %do%
#' @importFrom foreach %dopar%
# startingValuesFromModel = NULL
# starting.values.from.model = NULL
# startFromModel = NULL
# start.from.model

gpFitCV <- function (gp, fold.col, fold.fact = "1", folds = NULL, start.from.model = NULL,
	parallel = TRUE, fn.prefix = "", log.fn = if (parallel) "log-fold%f-%h-%p.txt" else NULL, dump.fn = if (parallel) "dump-fold%f-%h_%p"  else NULL,
	tr.max.lines = 5, pred.options = list(type = "terms", se.fit = TRUE),
	lmFit.options = list(lik.hyperpar = "fix"),
	lmFit.pred.options = list(type = "terms"),
	...)
{
	args <- list(fold.col = fold.col, fold.fact = fold.fact)

	if (is.null(gp[["data"]]))
		stop("The gp object doesn't contain scaled data; perpahs you need to call gpUnpack() on it (you can use compute = FALSE to save CPU time)")	
	
	# pre-check this, to provide more targetted error message if needed
	if (any(gp$hyperpar$component == ".lik") && (is.null(lmFit.options[["lik.hyperpar"]]) || lmFit.options[["lik.hyperpar"]] == "fix")) { 
		# there are likelihood hyperparameters, and to be fixed; check if lik.fix was specified
		if (is.null(lmFit.options[["lik.fix"]]))
			stop("missing the lik.fix argument in lmFit.options. It is needed since the likelihood template has hyperparameters, and you specified lmFit.options=list(lik.hyperpar = 'fix'). See ?lmFit for more details.")
	}

	fold.col <- gpFitCV__validate_and_get_fold_col(gp, fold.col, fold.fact)
		# validate fold.col and fold.fact arguments and get evaluated fold.col (as a vector)

	# now, init or check folds variable
	if (is.null(folds)) {
		folds <- sort(unique(fold.col))
	} else {
		stopifnot(all(folds %in% fold.col))
	}
	if (!is.null(start.from.model)) {
		stopifnot(class(start.from.model) == "gp")
		stopifnot(!is.null(start.from.model$fitCV))
		stopifnot(!is.null(start.from.model$fitCV$models))
		start.from.model.folds <- which(!sapply(start.from.model$fitCV$models, is.null))
			# we can't just do names(start.from.model$fitCV$models) since it is not named list, because the list indices are integers
		stopifnot(all(folds %in% start.from.model.folds))
	}
	
	ndigits <- floor(log10(max(fold.col))) + 1
	
	#if (fn.prefix != "") 
	#	fn.prefix <- paste0(
	
	if (parallel) {
		`%do_as_needed%` <- `%dopar%`
		if (!getDoParRegistered() || getDoParName() == "doSEQ")
			stop("for parallel = TRUE, parallel background for foreach() must be registered first")
	} else
		`%do_as_needed%` <- `%do%`
	
	wd <- getwd.keepsym()
	masterPID <- Sys.getpid()
	cat("1: getOption('warn') == ", getOption('warn'), "\n")
	cat(".packages(): \n")
	print(str(.packages()))
	#fold.run <- foreach (f = folds, .packages = c("gp")) %do_as_needed% {
	#fold.run <- foreach (f = folds, .packages = c("gp", "RTMB")) %do_as_needed% {	
	fold.run <- foreach (f = folds, .packages = .packages()) %do_as_needed% {
		if (getOption("warn") < 1)
			options(warn = 1) # tady veskere options() musi byt znova, protoze na worker/cluster se to neexportuje
		options(show.error.locations = TRUE)
		options(keep.source = TRUE)	
		options(traceback.max.lines = tr.max.lines)
		cat("2: getOption('warn') == ", getOption('warn'), "\n")
#library(gp)
#library(RTMB)		
		inv.logit <- RTMB::plogis # !!!! dirty hack!!! Shouldn't be needed!!!
		
		if (is.null(log.fn))
			log.fn2 <- NULL 
		else {
			log.fn2 <- gsub("%f", formatC(f, width = ndigits, flag = "0"), fixed = TRUE, log.fn)
			log.fn2 <- paste0(fn.prefix, log.fn2)
		}
		if (is.null(dump.fn))
			dump.fn2 <- NULL
		else {
			dump.fn2 <- gsub("%f", formatC(f, width = ndigits, flag = "0"), fixed = TRUE, dump.fn)
			dump.fn2 <- paste0(fn.prefix, dump.fn2)
		}
		parallelJobWrapper(working.dir = wd, masterPID = masterPID, log.fn = log.fn2, dump.fn = dump.fn2, parallel = parallel, tr.max.lines = tr.max.lines,
		{
			cat("3: getOption('warn') == ", getOption('warn'), "\n")
			gpcv <- gpPack(gp, maximum = TRUE)
			gpcv$fit <- NULL # delete the whole $fit object
			gpcv$fitCV <- NULL # don't forget also the $fitCV object!

			train_data <- gpDataSubset(gp$obsdata, fact = fold.fact, ind = (fold.col != f))
			test_data <- gpDataSubset(gp$obsdata, fact = fold.fact, ind = (fold.col == f))

			gpcv <- gp:::gpSetData(gpcv, train_data) # nechapu, proc tu najednou musi byt gp::: ... kdyz "gp" je mezi .packages! voser!

			if (!is.null(start.from.model))
				gpcv <- gpHyperparStartFromModel(gpcv, start.from.model$fitCV$models[[f]])

			# first fit & predict the null model - if something fails with an error, we want to fail it ASAP, not after some loong GP computation
			nullm <- do.call(lmFit, c(list(gp = gpcv, formula = ~1, data = gpcv$data), lmFit.options))
			test_data_lm <- gpDataPrepare(gpcv, test_data, scale.as = gpcv$data)
				# lmFit and predict.lmFit don't solve scaling, so we do it manually here, just in case...
				# but these scaling things might be worth checking in case some problems arise in models 
				# with complicated likelihood hyperparameters
			nullPredCV <- do.call(predict, c(list(nullm, newdata = test_data_lm), lmFit.pred.options))

			# now fit & predict the GP model - this might take a lot of time
			m <- gpFit(gpcv, ...)
			predCV <- do.call(predict, c(list(m, newdata = test_data), pred.options))
			
			m <- gpPack(m, maximum = TRUE)
			# pack it even more! These will be then restored in gpGetCVModel() :
			m$lik <- NULL
			m$obsdata <- NULL
			m$fit$a <- NULL
			if (!is.null(m$fit$stage1))
				m$fit$stage1$a <- NULL
			
			list(
				fold = f,
				m = m,
				predCV = predCV,
				nullPredCV = nullPredCV
			)
		})
	}
	args$folds <- folds
	gp$fitCV <- list(
		.pkg.version = pkg_build_info(),
		args = args,
		models = list(), # list of models for each fold
		predCV = NULL, # cross-validated prediction for the training dataset
		nullPredCV = NULL,
		stats = NULL # CV stats		
	)
	# put the results together
	# now, there might be two reasons to reindex the fold.col to the dimension of main table
	reindex.ff2main <- FALSE
	if (fold.fact != gp$GP_factor) { # first reason is this 
		# now, reindex the fold.col to the dimension of the prediction (gp$GP_factor)
		# thanks to the condition checked in gpFitCV__validate_and_get_fold_col(), the only case when reindexing might be needed is when GP_factor = "1" and fold.fact = something else (=something smaller)
		stopifnot(gp$GP_factor == "1") # consequence of the checks above
		stopifnot(fold.fact != "1") # consequence of the checks above
		# now, we have to reindex fold.col to the main table
		reindex.ff2main <- TRUE
	}
	if (gp$lik.reindex2main && gp$GP_factor != "1") { # this is the second reason
		# in this case, both fold.factor and GP_factor are something smaller than "1"
		# (the main table), and we want to reindex it to "1" = the main table
		stopifnot(fold.fact != "1")
		stopifnot(fold.fact == gp$GP_factor)
		reindex.ff2main <- TRUE
	}
	if (reindex.ff2main) {
		stopifnot(gpDataHasMainTable(gp[["obsdata"]])) # has to have it in this case
		fold_idx_col <- paste0(fold.fact, "_idx")
		fold.col <- fold.col[gp$obsdata[[1]][[fold_idx_col]]]
	}
	for (i in 1:length(fold.run)) { 
		f <- fold.run[[i]]$fold
		gp$fitCV$models[[f]] <- fold.run[[i]]$m
		
		# assemble GP predictions
		if (is.null(gp$fitCV$predCV))
			gp$fitCV$predCV <- as.data.frame(matrix(NA, nrow = length(fold.col), ncol = ncol(fold.run[[i]]$predCV), dimnames = list(NULL, colnames(fold.run[[i]]$predCV))))
		gp$fitCV$predCV[fold.col == f,] <- fold.run[[i]]$predCV

		# assemble null predictions
		if (is.null(gp$fitCV$nullPredCV))
			gp$fitCV$nullPredCV <- as.data.frame(matrix(NA, nrow = length(fold.col), ncol = ncol(fold.run[[i]]$nullPredCV), dimnames = list(NULL, colnames(fold.run[[i]]$nullPredCV))))
		gp$fitCV$nullPredCV[fold.col == f,] <- fold.run[[i]]$nullPredCV
		
	}
	
	# evaluate
	gp$fitCV[["stats"]] <- tryNull(gp$lik$predPerf(data = gp$data, pred = gp$fitCV$predCV, pred.null = gp$fitCV$nullPredCV))
	
	gp
}

# internal helper function that validates fold.col and fold.fact arguments and returns evaluated fold.col (as a vector)
gpFitCV__validate_and_get_fold_col <- function (gp, fold.col, fold.fact)
{
	if (gp$GP_factor != "1")
		if (fold.fact != gp$GP_factor)
			stop("If the factor corresponding to the GP (gp$GP_factor) is some real grouping factor (!= \"1\"), then the fold.fact must be that factor")
	if (is.character(fold.col) && length(fold.col) == 1) { # it is a name of a column in the main table
		if (fold.fact != "1")
			stop("Specifying fold.col by column name only works for fold.fact == '1'")
		stopifnot(gpDataHasMainTable(gp$obsdata))
		if (!fold.col %in% colnames(gp$obsdata[[1]]))
			stop("column ", fold.col, " does not exist in the main table of the gp$obsdata.")
		fold.col <- gp$obsdata[[1]][[fold.col]]
	}
	# fold.col is now a vector along fold.fact
	# check it now
	stopifnot(length(fold.col) == gpDataSize(gp$obsdata, fold.fact))
	# it should be a vector of integer numbers from 1 to N, where N is the number of folds
	stopifnot(is.integer(fold.col))
	stopifnot(all(1:max(fold.col) == sort(unique(fold.col))))
	return(fold.col)
}

#' Get the cross-validation model for given fold
#'
#' @param gp GP model object
#' @param fold integer, number of the cross-validation fold
#' @return GP model object for the fold \code{fold}
#' @export
gpGetCVModel <- function(gp, fold)
{
	stopifnot(!is.null(gp[["fitCV"]]))
	stopifnot(!is.null(gp$fitCV[["models"]]))
	stopifnot(!is.null(gp$fitCV[["args"]]))
	stopifnot(is.numeric(fold))
	stopifnot(fold >= 1)
	stopifnot(fold <= length(gp$fitCV$models))
	
	fold.fact <- gp$fitCV$args$fold.fact
	fold.col <- gpFitCV__validate_and_get_fold_col(gp, gp$fitCV$args$fold.col, gp$fitCV$args$fold.fact)
		# validate fold.col and fold.fact arguments and get evaluated fold.col (as a vector)
	
	m <- gp$fitCV$models[[fold]]
	m$lik <- gp$lik
	m$obsdata <- gpDataSubset(gp$obsdata, fact = fold.fact, ind = (fold.col != fold))
	
	m <- gpUnpack(m, compute = FALSE)
	# would need to call gpFitLaplace() or something to calculate `a`, we don't even have that!
	
	return(m)
}

gpGetFoldCol <- function(gp)
{
	stopifnot(!is.null(gp[["fitCV"]]))
	stopifnot(!is.null(gp$fitCV[["args"]]))
	fold.fact <- gp$fitCV$args$fold.fact
	fold.col <- gpFitCV__validate_and_get_fold_col(gp, gp$fitCV$args$fold.col, gp$fitCV$args$fold.fact)
		# validate fold.col and fold.fact arguments and get evaluated fold.col (as a vector)
	return (fold.col)
}

