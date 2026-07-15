#' Perform an N-fold cross-validation of a linear model fitted on a GP model object
#'
#' Function will fit the linear model for each of the N cross-validation folds, and evaluate cross-validation prediction for the whole training data set.
#'
#' @param gp GP model object
#' @param fold.col either a name of a column in the main table, or a vector along factor \code{fold.fact}.
#'		The column (or the supplied vector) must be a vector of integers from \code{1} to \code{N} (\code{N} being the number of folds), 
#'		specifying the number of the cross-validation fold the given record belongs to.
#' @param fold.fact a factor along which the cross-validation folds (\code{fold.col}) are specified. Factor \code{"1"} means the folds are specified for the rows of the main table.
#'		Note that if the GP dimension is given by some real grouping factor (i.e. \code{gp$GP_factor != "1"}), then the \code{fold.fact} must be that factor.
#' @param folds integer vector of folds to fit; if \code{NULL} (the default), all folds are fit
#' @param parallel should the cross-validation models be run in parallel? Uses foreach(), and requires parallel background for foreach() to be already registered.
#' @param fn.prefix character; prefix of file names for log files (\code{log.fn}) and dump files (\code{dump.fn}), in case these are not \code{NULL}
#' @param log.fn if not \code{NULL}, standard and error output of each fold model job will be saved into a log file with this file name. Special \code{\%} sequences can be used, see Details below. 
#' @param dump.fn if not \code{NULL}, debug dump of given fold model will be saved upon an error, with this file name (the .rda extension will be added to it). Special \code{\%} sequences can be used, see Details below. 
#' @param tr.max.lines the \code{max.lines} parameter for \code{traceback}, i.e. the maximum number of lines printed per call when error occurs
#' @param pred.options list of options to be passed to the \code{predict.lmFit} method.
#' @param ... options to be passed to the \code{\link{lmFit}} method.
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

lmFitCV <- function (gp, fold.col, fold.fact = "1", folds = NULL,
	parallel = FALSE, fn.prefix = "", log.fn = NULL, dump.fn = NULL,
	tr.max.lines = 5, pred.options = list(type = "terms", se.fit = TRUE),
	...)
{
	args <- list(fold.col = fold.col, fold.fact = fold.fact)
	lm <- lmFit(gp, ...)
	lmFit.args <- list(...)
	
# note that we don't use the data from the gp object any more; we use the data
# provided by the data argument in the lmFitCV() call

	fold.col <- gpFitCV__validate_and_get_fold_col(gp, fold.col, fold.fact, data = lm$data)
		# validate fold.col and fold.fact arguments and get evaluated fold.col (as a vector)

	# now, init or check folds variable
	if (is.null(folds)) {
		folds <- sort(unique(fold.col))
	} else {
		stopifnot(all(folds %in% fold.col))
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
	#cat("1: getOption('warn') == ", getOption('warn'), "\n")
	#cat(".packages(): \n")
	#print(str(.packages()))
	#fold.run <- foreach (f = folds, .packages = c("gp")) %do_as_needed% {
	#fold.run <- foreach (f = folds, .packages = c("gp", "RTMB")) %do_as_needed% {	
	fold.run <- foreach (f = folds, .packages = .packages()) %do_as_needed% {
		if (getOption("warn") < 1)
			options(warn = 1) # tady veskere options() musi byt znova, protoze na worker/cluster se to neexportuje
		options(show.error.locations = TRUE)
		options(keep.source = TRUE)	
		options(traceback.max.lines = tr.max.lines)
		#cat("2: getOption('warn') == ", getOption('warn'), "\n")
		
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
			#cat("3: getOption('warn') == ", getOption('warn'), "\n")

			train_data <- gpDataSubset(lm$data, fact = fold.fact, ind = (fold.col != f))
			test_data <- gpDataSubset(lm$data, fact = fold.fact, ind = (fold.col == f))
				# !!! lmFit and predict.lmFit don't solve scaling
				# but these scaling things might be worth checking in case some problems arise in models 
				# with complicated likelihood hyperparameters

			# fit & predict the null model
			nullFit.args <- lmFit.args
			nullFit.args$gp <- NULL
			nullFit.args$formula <- NULL
			nullFit.args$data <- NULL
			nullm <- do.call(function(...) lmFit(gp, formula = ~1, data = train_data, ...), nullFit.args)
			nullPredCV <- do.call(function(...) predict(nullm, newdata = test_data, ...), pred.options)

			# now fit & predict the linear model
			fit.args <- lmFit.args
			fit.args$gp <- NULL
			fit.args$data <- NULL
			m <- do.call(function(...) lmFit(gp, data = train_data, ...), fit.args)
			predCV <- do.call(function(...) predict(m, newdata = test_data, ...), pred.options)
			
			list(
				fold = f,
				m = m,
				predCV = predCV,
				nullPredCV = nullPredCV
			)
		})
	}
	args$folds <- folds
	lm$fitCV <- list(
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
		stopifnot(gpDataHasMainTable(lm[["data"]])) # has to have it in this case
		fold_idx_col <- paste0(fold.fact, "_idx")
		fold.col <- fold.col[lm$data[[1]][[fold_idx_col]]]
	}
	for (i in 1:length(fold.run)) { 
		f <- fold.run[[i]]$fold
		#lm$fitCV$models[[f]] <- fold.run[[i]]$m # do not store the models in the lmFitCV() case, can be too big for no reason
		
		# assemble linear model predictions
		if (is.null(lm$fitCV$predCV))
			lm$fitCV$predCV <- as.data.frame(matrix(NA, nrow = length(fold.col), ncol = ncol(fold.run[[i]]$predCV), dimnames = list(NULL, colnames(fold.run[[i]]$predCV))))
		lm$fitCV$predCV[fold.col == f,] <- fold.run[[i]]$predCV

		# assemble null predictions
		if (is.null(lm$fitCV$nullPredCV))
			lm$fitCV$nullPredCV <- as.data.frame(matrix(NA, nrow = length(fold.col), ncol = ncol(fold.run[[i]]$nullPredCV), dimnames = list(NULL, colnames(fold.run[[i]]$nullPredCV))))
		lm$fitCV$nullPredCV[fold.col == f,] <- fold.run[[i]]$nullPredCV
		
	}
	
	# evaluate
	lm$fitCV[["stats"]] <- tryNull(gp$lik$predPerf(data = lm$data, pred = lm$fitCV$predCV, pred.null = lm$fitCV$nullPredCV))
	
	lm
}
