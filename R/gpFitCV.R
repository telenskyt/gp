#' Perform an N-fold cross-validation of a Gaussian process model
#'
#' Function will fit the GP model for each of the N cross-validation folds, and evaluate cross-validation prediction for the whole training data set.
#'
#' @param gp GP model object
#' @param fold.col either a name of a column in the main table, or a vector along factor \code{fold.fact}.
#'		The column (or the supplied vector) must be a vector of integers from \code{1} to \code{N} (\code{N} being the number of folds), 
#'		specifying the number of the cross-validation fold the given record belongs to.
#' @param fold.fact a factor along which the cross-validation folds (\code{fold.col}) are specified. Factor \code{"1"} means the folds are specified for the rows of the main table.
#'		Note that if the GP dimension is given by some real grouping factor (i.e. \code{p$GP_factor != "1"}), then the \code{fold.fact} must be that factor.
#' @param folds integer vector of folds to fit; if \code{NULL} (the default), all folds are fit
#' @param start.from.model an object of class \code{gp} to take the starting values from. If specified, then for each cross-validation fold model to be fit, 
#' 		the starting values of the hyperparameters will be taken from the corresponding fold model wherever possible.
#' @param parallel should the cross-validation models be run in parallel? Uses foreach(), and requires parallel background for foreach() to be already registered.
#' @param fn.prefix character; prefix of file names for log files (\code{log.fn}) and dump files (\code{dump.fn}), in case these are not \code{NULL}
#' @param log.fn if not \code{NULL}, standard and error output of each fold model job will be saved into a log file with this file name. Special \code{\%} sequences can be used, see Details below. 
#' @param dump.fn if not \code{NULL}, debug dump of given fold model will be saved upon an error, with this file name (the .rda extension will be added to it). Special \code{\%} sequences can be used, see Details below. 
#' @param ... options to be passed to the \code{\link[gpFit]{gpFit()}} method.
#'
#' @details The arguments \code{log.fn} and \code{dump.fn} allow for special sequences:
#' - \code{\%f} - number of the cross-validation fold
#' - \code{\%h} - hostname, i.e. the name of the machine where the worker job runs
#' - \code{\%p} - process ID of the worker job
#'
#' @export
# startingValuesFromModel = NULL
# starting.values.from.model = NULL
# startFromModel = NULL
# start.from.model

gpFitCV <- function (gp, fold.col, fold.fact = "1", folds = NULL, start.from.model = NULL,
	parallel = TRUE, fn.prefix = "", log.fn = "log-fold%f-%h-%p.txt", dump.fn = "dump-fold%f-%h_%p",
	...)
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
	#fold.run <- foreach (f = folds, .packages = c("gp")) %do_as_needed% {
	#fold.run <- foreach (f = folds, .packages = c("gp", "RTMB")) %do_as_needed% {	
	fold.run <- foreach (f = folds, .packages = .packages()) %do_as_needed% {	
		options(warn = 2) # tady veskere options() musi byt znova, protoze na worker/cluster se to neexportuje
		options(show.error.locations = TRUE)
		options(keep.source = TRUE)	
		
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
		parallelJobWrapper(working.dir = wd, masterPID = masterPID, log.fn = log.fn2, dump.fn = dump.fn2, 
		{
			gpcv <- gpPack(gp, maximum = TRUE)
			gpcv$fit <- NULL # delete the whole $fit object

			train_data <- gpDataSubset(gp$obsdata, fact = fold.fact, ind = (fold.col != f))
			test_data <- gpDataSubset(gp$obsdata, fact = fold.fact, ind = (fold.col == f))

			{ ### Integrate the newly subsetted data into the CV model
			  ### !! This code is duplicated with gp() - perhaps a sign that the interface should have been different, not passing data to gp() but to gpFit() and gpFitCV()?
			gpcv$obsdata <- train_data
			gpcv$data <- gpDataPrepare(gpcv, gpcv$obsdata)

			gpcv$covComp_df <- gp:::gpComponentsTable(gpcv)

			gp_size <- gpDetermineSize(gpcv)
			gpcv$GP_size <- gp_size$size
			stopifnot(gpcv$GP_factor == gp_size$fact) # bude vzdy character string ruzny od "", NA, NULL, viz gpDetermineSize()
			}
			
			if (!is.null(start.from.model))
				gpcv <- gpHyperparStartFromModel(gpcv, start.from.model$fitCV$models[[f]])
					
			m <- gpFit(gpcv, ...)
			predCV <- predict(m, test_data, type = "latent", se.fit = TRUE)
			
			m <- gpPack(m, maximum = TRUE)
			# pack it even more! :
			m$obsdata <- NULL
			m$fit$a <- NULL
			
			list(
				fold = f,
				m = m,
				predCV = predCV
			)
		})
	}
	gp$fitCV <- list(
		models = list(), # list of models for each fold
		predCV = NULL, # cross-validated prediction for the training dataset
		stats = NULL # CV stats		
	)
	# put the results together
	# now, reindex the fold.col to the dimension of the prediction (gp$GP_factor)
	# thanks to the condition above, the only case when reindexing might be needed is when GP_factor = "1" and fold.fact = something else.
	if (fold.fact != gp$GP_factor) {
		stopifnot(gp$GP_factor == "1") # consequence of the checks above
		stopifnot(fold.fact != "1") # consequence of the checks above
		# now, we have to reindex fold.col to the main table
		stopifnot(gpDataHasMainTable(gp$obsdata)) # has to have it in this case
		fold_idx_col <- paste0(fold.fact, "_idx")
		fold.col <- fold.col[gp$obsdata[[1]][[fold_idx_col]]]
	}
	for (i in 1:length(fold.run)) { 
		f <- fold.run[[i]]$fold
		gp$fitCV$models[[f]] <- fold.run[[i]]$m
		
		if (is.null(gp$fitCV$predCV))
			gp$fitCV$predCV <- as.data.frame(matrix(NA, nrow = length(fold.col), ncol = ncol(fold.run[[i]]$predCV), dimnames = list(NULL, colnames(fold.run[[i]]$predCV))))
		gp$fitCV$predCV[fold.col == f,] <- fold.run[[i]]$predCV
		
	}
	#gp$fitCV$stats <- ... !!!
	gp
}

