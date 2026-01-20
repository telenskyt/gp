

#' Redirect stdout and stderr to a log file
#'
#' Redirects R's standard output and messages (stderr)
#' to the specified file. On Unix systems, can optionally keep output visible on the
#' console using \code{tee} and add timestamps via the \code{ts} command.
#'
#' @param filename character. Path to the log file.
#' @param tee logical. If \code{TRUE}, stream output to console as well (Unix only). currently on windows, STDERR cannot be split (i.e. will not be displayed on the console).
#' @param all_or_nothing logical. Makes only sense for tee = TRUE. When TRUE, perform "tee" only on systems where full funcionality is available 
#					(on Windows do nothing now)
# 				- when FALSE, do tee as much as you can (on Windows, STDERR won't be displayed)
#' @param add_timestamp logical. Prepend timestamps using \code{ts} (Unix only, requires `ts` from moreutils).
#' @param append logical. If \code{TRUE}, append to \code{filename} instead of overwriting.
#' @param autoWarn1 logical. If \code{TRUE}, it will set \code{options(warn = 1)} so that warnings are printed to the log immediately as well.
#'
#' @return A connection object used for logging, or \code{NULL} on failure.
#'
#' @examples
#' \dontrun{
#' LOG <- startLog("mylog.txt", tee = interactive(), append = TRUE)
#' # ... do work ...
#' endLog(LOG)
#' }
#'
#' @export

# my orig help:
# redirect STDOUT and STDERR to log file `filename`
# tee - when TRUE,  keep it flowing to the console (tee). This now only works on unix systems,
# 		currently on windows STDERR cannot be split (i.e. will not be displayed)
# all_or_nothing - makes only sense for tee = TRUE. When TRUE, perform "tee" only on systems where full funcionality is available 
#					(on Windows do nothing now)
# 				- when FALSE, do tee as much as you can (on Windows, STDERR won't be displayed)
# add_timestamp - will work only on UNIX, and requires `ts` from moreutils
# autoWarn1 - when user is logging things, he probably assumes that warnings will be there as well. So we make sure they are (autoWarn1 = TRUE).
startLog <- function (filename, tee = FALSE, all_or_nothing = TRUE, add_timestamp = FALSE, append = FALSE, autoWarn1 = TRUE)
{
	if (append)
		mode <- "at"
	else
		mode <- "wt"
	conn1 <- NULL
	if (autoWarn1 && getOption("warn") == 0) {
		options(warn = 1)
		#warning("options(warn) was 0, setting to 1")
	}
	if (!tee) {
		if (add_timestamp && .Platform$OS.type == "unix") {
			if (append)
				conn1 <- pipe(paste0("ts >> ", filename), open = "wt")
			else
				conn1 <- pipe(paste0("ts > ", filename), open = "wt")
		} else
			conn1 <- file(filename, open = mode)
		sink(conn1)
		sink(conn1, type = "message")
# from now on, tee = TRUE
	} else if (.Platform$OS.type == "unix") { # pro veskery non-parallel veci; vse chci logovat; nevim jak jinak nez vse dohromady
										# warnings, errors, '!!!' etc si z toho vygrepuju sam
		if (add_timestamp) {
			if (append)
				conn1 <- pipe(paste0("ts | tee -a ", filename), open = "wt")
			else
				conn1 <- pipe(paste0("ts | tee ", filename), open = "wt")
		} else {
			if (append)
				conn1 <- pipe(paste0("tee -a ", filename), open = "wt")
			else
				conn1 <- pipe(paste0("tee ", filename), open = "wt")
		}
		sink(conn1)
		sink(conn1, type = "message")
		# Q1: will this be passed on the child parallel things?
	} else if (!all_or_nothing) { 
		conn1 <- file(filename, open = mode)
		sink(conn1, split = TRUE)
		sink(conn1, type = "message") # nejsou split, nebudou videt na konzoli!!!
	}
	conn1
}

#
# conn - connection returned by startLog()

#' Stop logging started by startLog()
#'
#' Stop redirecting the stdout and stderr to a log file and close the connection returned by
#' \code{startLog}.
#'
#' @param conn Connection object returned by \code{startLog()}.
#' 
#' @examples
#' \dontrun{
#' LOG <- startLog("mylog.txt", tee = interactive(), append = TRUE)
#' # ... do work ...
#' endLog(LOG)
#' }
#' @export
endLog <- function (conn)
{
	if (!is.null(conn)) {
		sink(type = "message") 
		sink()
		close(conn)
	}
}


#' Get working directory as user specified it (keeping symbolic links in the path)
#'
#' Like getwd(), but keeps the path originally specified by the user, with symbolic links intact. This is a difference from \code{getwd()}, because on UNIX systems the path 
#' might contain symlinks, and getwd() will then return a different path, with those symlinks resolved.
#' @export
getwd.keepsym <- function()
{
	if (.Platform$OS.type == "unix") {
		wd <- Sys.getenv('PWD')
		stopifnot(wd != "")
	} else {
		wd <- getwd()
	}
	wd
}




expand_special <- function (fn)
{
	hostname <- Sys.info()["nodename"]
	#hostname.abbrev <- sub("telensky-vypocty-?", "", hostname)
	#if (hostname.abbrev == "") hostname.abbrev <- "1"	
	#fn <- gsub("%H", hostname.abbrev, fixed = TRUE, fn)
	fn <- gsub("%h", hostname, fixed = TRUE, fn)
	fn <- gsub("%p", Sys.getpid(), fixed = TRUE, fn)
	fn
}



#' Wrapper for parallel/background jobs
#'
#' Useful wrapper for parallel jobs, run by e.g. foreach(), or, in general, for any non-interactive (background) jobs. 
#' It will make sure the working directory is set correctly, standard and error outputs are logged, and in case of 
#' error, a dump file will be saved for debugging purposes.
#'
#' @param expr expression - the job code to be run
#' @param working.dir working directory
#' @param masterPID the PID (process ID) of the process that is dispatching this parallel computation
#' @param log.fn if not \code{NULL}, standard and error output of the job will be saved into a log file with this file name. Special \code{\%} sequences can be used, e.g. \code{"log-\%h_\%p.txt"}, see Details below. 
#' @param dump.fn if not \code{NULL}, debug dump will be saved upon an error, with this file name (the .rda extension will be added to it). Special \code{\%} sequences can be used, e.g. \code{"dump-\%h_\%p"}, see Details below. 
#' @param parallel should be set to \code{TRUE} for the wrapper to work normally. The \code{FALSE} setting is just for some perhaps
#' debugging purposes, when we want to run the job in non-parallel, interactive mode, in which case most of the wrapper function (logging, setting working directory)
#' will be disabled.
#' @details The arguments \code{log.fn} and \code{dump.fn} allow for special sequences:
#' - \code{\%h} - hostname, i.e. the name of the machine where the worker job runs
#' - \code{\%p} - process ID of the worker job
#' @export
parallelJobWrapper <- function (expr, working.dir = NULL, masterPID = NULL, log.fn = NULL, dump.fn = NULL, parallel = TRUE)

#myParallel
#jobWrapper
#bgJobWrapper
#bgJob
#parallelJob
#parallelJobWrapper
#parJobWrapper

{
	res <- withCallingHandlers({
		if (parallel) {
			hostname <- Sys.info()["nodename"]
			if (!is.null(working.dir)) {
				setwd(working.dir)
			}						
			# spust logovani - dej logu do jmena hostname a PID, kdyby nahodou na tom samem tasku makalo vic workeru (u doRedis clovek nikdy nevi!)
			if (!is.null(log.fn)) {
				log.fn <- expand_special(log.fn)
				fold.LOG <- startLog(filename = log.fn, add_timestamp = TRUE)
				on.exit({ cat("closing log file\n"); endLog(fold.LOG); cat("(after closing log, shouldn't be seen in any log)\n"); }) 
					# we have to do it like this, because of the doRedis bug (worker will not close the sink; po case to zpusobi chybu "sink stack is full". Not reported yet)
			}
			cat("Running in parallel: my PID = ", Sys.getpid(), " at ", hostname, "; masterPID = ", masterPID, "\n")
			cat("Starting at the following memory: gc(reset = TRUE):\n")
			message("testing stderr (this goes to message connection, which is like stderr)\n")
			print(gc(reset = TRUE))
			cat("\n\n")
		}
		
		expr
	}, error = function (e) {
		#traceback() # i zde vyhodi no traceback available
		cat("error: ")
		print(e)
		#cat("\n")
		
		if (parallel) {
			cat("The stack:\n") # nevim jak ho vypsat, nic nefacha:
			#recover() # nic nevypise
			#browser() # nic nevypise
			#traceback() # vyhodi no traceback available! # trochu me stve, ze to nedava ten samy stack format jako recover()... mozna todo: vykuchat to z recover()...
			traceback(sys.calls()) # this works!!! 
		}
		
		if (!is.null(dump.fn)) {
			dump.fn <- expand_special(dump.fn)
			cat("dumping stack, variables etc. to ", dump.fn, "\n")
			dump.frames(dumpto = dump.fn, to.file = TRUE, include.GlobalEnv = TRUE) # use debugger(loadVar("*", dump.fn)) to debug...
		}
	})
	gc() # free up the memory because e.g. doRedis workers don't call gc() after it finishes (they should!!)
	return(res)
}

#' Enhanced version of foreach()
#' 
#' Enhanced version of \code{foreach()} that takes care of setting correct working directory, logging standard and error output to a log file, and dumping debug info upon error.
#' 
#' @param .pass.wd should the worker process set the same working directory as is on the master, before the job is executed? In most cases, \code{TRUE} will be the desired value (default).
#'		Keep in mind that in some setups, e.g. with \code{library(doRedis)}, one might have the workers running on different machine.
#' @param .working.dir working directory to be set before the job is executed. Will override \code{.pass.wd} if set.
#' @param .log.fn if not \code{NULL}, standard and error output of the job will be saved into a log file with this file name. 
#' 		Special \code{\%} sequences can be used, e.g. \code{"log-\%h_\%p.txt"}, see Details below. 
#'		Note that this parameter is evaluated in the context of the inner loop, 
#' (!!! o to jsem se snazil ale to nefunguje bohuzel), so you might happily use the 
#' 		foreach2() iterator variables in the expression!
#' @param .dump.fn if not \code{NULL}, debug dump will be saved upon an error, with this file name (the .rda extension will be added to it). 
#' 		Special \code{\%} sequences can be used, e.g. \code{"dump-\%h_\%p"}, see Details below.
#' Note that this parameter is evaluated in the context of the inner loop (!!! o to jsem se snazil ale to nefunguje bohuzel), 
#'	so you might happily might use the \code{foreach2()} iterator variables in the expression!
#' @param ... arguments passed to \code{\link[foreach]{foreach}()}.
#' @details The arguments \code{log.fn} and \code{dump.fn} allow for special sequences:
#' - \code{\%h} - hostname, i.e. the name of the machine where the worker job runs
#' - \code{\%p} - process ID of the worker job
#' @export
foreach2 <- function (.pass.wd = TRUE, .working.dir = NULL, .log.fn = NULL, .dump.fn = NULL, .parallel = TRUE, ...)
{
	obj <- foreach(...)
	masterPID <- Sys.getpid()
	if (is.null(.working.dir) && .pass.wd) {
		.working.dir <- getwd.keepsym()
	}
	obj$foreach2 <- list(
		parallel = .parallel,
		working.dir = .working.dir,
		masterPID = masterPID,
		log.fn = .log.fn, # evaluate this in the context of the loop
		dump.fn = .dump.fn # evaluate this in the context of the loop
		#log.fn = quote(.log.fn), # evaluate this in the context of the loop
		#dump.fn = quote(.dump.fn) # evaluate this in the context of the loop
	)
	obj
}

#' Enhanced version of `%dopar%`
#' 
#' 	Enhanced version of \code{\%dopar\%} that works with \code{foreach2()} and takes care of setting correct working directory, logging standard and error output to a log file, and dumping debug info upon error.
#' @export
#
#
# tato implementace sice funguje, ale jen bohuzel prasacky, ze jsem zkopiroval a upravil kod z foreach::`%dopar%`... a navic by to nefungovalo pro `%do_as_needed%` <- `%do%`
`%dopar2%` <- function (obj, ex)
{
	args <- obj$foreach2
	if (args$parallel) 
		`%do_as_needed%` <- `%dopar%`
	else
		`%do_as_needed%` <- `%do%`

	if (0) {
	obj %do_as_needed% {
		args <- c(args, list(expr = quote(ex)))
		do.call(parallelJobWrapper, args)
	}
	}
	
	#obj %do_as_needed% quote(ex) # doesn't work either
	#obj %do_as_needed% substitute(ex, parent.frame())
	#obj %dopar% ex

    e <- foreach:::getDoPar()
    e$fun(obj, substitute({
		#ex
		#args <- c(args, list())
		gp::parallelJobWrapper(
			parallel = args$parallel,
			working.dir = args$working.dir,
			masterPID = args$masterPID,
			log.fn = args$log.fn, # evaluate this in the context of the loop
			dump.fn = args$dump.fn, # evaluate this in the context of the loop		
			expr = ex)
	}), parent.frame(), e$data)	
	#e$fun(obj, substitute(ex), parent.frame(), e$data)	
}
