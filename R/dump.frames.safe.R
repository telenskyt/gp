
# Internal variant of utils::dump.frames() that avoids passing raw calls to
# limitedLabels(), because that can fully deparse huge do.call()-generated calls.
# (This is a workaround for R problem, should be reported as a bug)
#
# This intentionally avoids using utils:::limitedLabels() as an internal helper.
# Calls are first deparsed with a bounded number of lines, then labels are
# trimmed locally using the same width limits used by limitedLabels().
#
# Created using ChatGPT & Codex

.dump.frames.safe <- function(dumpto = "last.dump", to.file = FALSE, include.GlobalEnv = FALSE,
	max.lines = getOption("gp.dump.frames.max.lines", 1L),
	compress = getOption("gp.dump.frames.compress", FALSE))
{
	max.lines <- suppressWarnings(as.integer(max.lines)[1L])
	if (is.na(max.lines) || max.lines < 1L)
		max.lines <- 1L

	calls <- sys.calls()
	last.dump <- sys.frames()

	safe.labels <- vapply(calls, function(cl) {
		paste(deparse(cl, nlines = max.lines), collapse = " ")
	}, character(1L))

	maxwidth <- getOption("width") - 5L
	if (is.null(maxwidth) || maxwidth < 40L)
		maxwidth <- 40L
	maxwidth <- min(maxwidth, 1000L)

	names(last.dump) <- strtrim(safe.labels, maxwidth)

	if (include.GlobalEnv) {
		last.dump <- c(
			".GlobalEnv" = as.environment(as.list(.GlobalEnv, all.names = TRUE)),
			last.dump
		)
	}

	last.dump <- last.dump[-length(last.dump)]
	attr(last.dump, "error.message") <- geterrmessage()
	class(last.dump) <- "dump.frames"

	if (dumpto != "last.dump")
		assign(dumpto, last.dump)

	if (to.file) {
		save(
			list = dumpto,
			file = paste0(dumpto, ".rda"),
			envir = environment(),
			compress = compress
		)
	} else {
		assign(dumpto, last.dump, envir = .GlobalEnv)
	}

	invisible()
}
