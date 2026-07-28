# Return the optimization stages stored in a fitted gp object.
.gpOptimizationStages <- function(x)
{
	fit <- x[["fit"]]
	if (is.null(fit))
		stop("Model object has not been fit yet: you need to call gpFit() first")

	if (!is.null(fit[["stage1"]]))
		return(list("Stage 1" = fit$stage1, "Stage 2" = fit))

	stage <- "Optimization"
	if (isTRUE(fit$args$two.stage)) {
		if (2 %in% fit$args$stages)
			stage <- "Stage 2"
		else if (1 %in% fit$args$stages)
			stage <- "Stage 1"
	}
	setNames(list(fit), stage)
}


# Prepare one stage of the iteration table.
.gpIterTable <- function(fit, use.prior, factr)
{
	info <- fit[["iter_info"]]
	if (is.null(info))
		return(NULL)

	if (use.prior) {
		objective <- info$negLogPosterior
		objective.name <- "negLogPosterior"
	} else {
		objective <- info$nlml
		objective.name <- "nlml"
	}

	objective.prev <- c(NA, head(objective, -1))
	obj.diff <- objective - objective.prev
	obj.rel.chg <- obj.diff / pmax(objective, objective.prev, 1)

	factr2 <- factr
	show.factr <- rep("", length(obj.rel.chg))
	while (factr2 <= 1e15) {
		cond <- -obj.rel.chg <= .Machine$double.eps * factr2 &
			-obj.rel.chg >= 0
		cond[is.na(cond)] <- FALSE
		if (any(cond)) {
			first <- which(cond)[1]
			if (show.factr[first] == "")
				show.factr[first] <- paste0("< dbl.eps * 1e",
					log(factr2) / log(10))
			if (first <= 2)
				break
		}
		factr2 <- factr2 * 10
	}

	keep <- c("iter", "LA_iterations", "nlml", "negLogPosterior",
		"lastLAObjDiff", "fStartReset", "AUC", "time")
	out <- info[, keep, drop = FALSE]

	returned <- rep(FALSE, nrow(out))
	if (!is.null(fit[["iter_h"]]) &&
		nrow(fit$iter_h) == nrow(out) &&
		length(fit[["h"]]) == ncol(fit$iter_h)) {
		returned <- apply(fit$iter_h, 1, function(h)
			isTRUE(all.equal(as.numeric(h), as.numeric(fit$h),
				check.attributes = FALSE)))
		if (any(returned))
			returned[-max(which(returned))] <- FALSE
	}
	out$iter <- as.character(out$iter)
	out$iter[returned] <- paste0(out$iter[returned], "*")

	out$obj.diff <- obj.diff
	out$obj.rel.chg <- obj.rel.chg
	out$factr <- show.factr
	attr(out, "objective") <- objective.name
	attr(out, "returned") <- any(returned)
	out
}


#' Display hyperparameter optimization iterations
#'
#' Displays the recorded objective evaluations and Laplace-approximation
#' information for each optimization stage. An asterisk after \code{iter}
#' marks the evaluation whose hyperparameters were returned by the optimizer.
#'
#' @param x A fitted object of class \code{gp}.
#' @param ... Additional arguments, currently ignored.
#'
#' @return Invisibly returns the displayed data frame. For a two-stage fit,
#' the returned data frame also contains a \code{stage} column.
#'
#' @export
iter <- function(x, ...)
{
	UseMethod("iter")
}


#' @rdname iter
#' @method iter gp
#' @export
iter.gp <- function(x, ...)
{
	stages <- .gpOptimizationStages(x)
	fit <- x[["fit"]]

	use.prior <- TRUE
	if (!is.null(fit$args$use.prior))
		use.prior <- fit$args$use.prior

	factr <- 1e7
	if (!is.null(fit$args$opt.control$factr))
		factr <- fit$args$opt.control$factr

	result <- list()
	returned.marked <- FALSE
	for (stage in names(stages)) {
		tab <- .gpIterTable(stages[[stage]], use.prior, factr)
		if (is.null(tab))
			next

		if (length(stages) > 1 || stage != "Optimization")
			cat(stage, "\n")
		print(tab)
		cat("\t\t(objective =", attr(tab, "objective"), ")\n\n")

		returned.marked <- returned.marked || attr(tab, "returned")
		tab$stage <- stage
		tab <- tab[, c("stage", setdiff(names(tab), "stage")), drop = FALSE]
		result[[stage]] <- tab
	}

	if (!length(result))
		stop("The fitted model does not contain hyperparameter optimization iterations.")

	if (returned.marked)
		cat("* hyperparameters returned by the optimizer\n")
	cat("factr =", factr, ", factr * double.eps =",
		factr * .Machine$double.eps, "\n")

	out <- do.call(rbind, result)
	rownames(out) <- NULL
	invisible(out)
}


#' Check model convergence
#'
#' Gives a short convergence verdict based on the optimizer status and the
#' stored Laplace-approximation iterations.
#'
#' @inheritParams iter
#'
#' @return Invisibly returns \code{TRUE} when all stored optimization stages
#' and Laplace approximations converged, and \code{FALSE} otherwise.
#'
#' @export
convergence <- function(x, ...)
{
	UseMethod("convergence")
}


#' @rdname convergence
#' @method convergence gp
#' @export
convergence.gp <- function(x, ...)
{
	stages <- .gpOptimizationStages(x)
	optimizer.lines <- character()
	optimizer.ok <- TRUE
	optimized <- FALSE
	LA.at.max <- 0
	LA.increased <- 0

	for (stage in names(stages)) {
		fit <- stages[[stage]]

		if (!is.null(fit[["optim"]])) {
			optimized <- TRUE
			code <- fit$optim$convergence
			ok <- identical(code, 0L) || identical(code, 0)
			optimizer.ok <- optimizer.ok && ok
			label <- if (length(stages) > 1 || stage != "Optimization")
				paste0(stage, " optimizer")
			else
				"Hyperparameter optimizer"
			message <- fit$optim$message
			if (is.null(message) || !nzchar(message))
				message <- "no convergence message"
			optimizer.lines <- c(optimizer.lines,
				paste0(label, ": ", if (ok) "OK" else "NOT OK",
					" (code ", code, ") - ", message))
		}

		if (!is.null(fit[["iter_info"]])) {
			tol <- if (is.null(fit[["tol"]])) 1e-6 else fit$tol
			stage.at.max <- 0
			stage.increased <- 0
			if (!is.null(fit[["itmax"]]))
				stage.at.max <-
					sum(fit$iter_info$LA_iterations >= fit$itmax, na.rm = TRUE)
			if ("lastLAObjDiff" %in% names(fit$iter_info))
				stage.increased <-
					sum(fit$iter_info$lastLAObjDiff > tol, na.rm = TRUE)
		} else {
			tol <- if (is.null(fit[["tol"]])) 1e-6 else fit$tol
			stage.at.max <- 0
			stage.increased <- 0
		}

		if (!is.null(fit[["itmax"]]) && !is.null(fit[["iterations"]]) &&
			fit$iterations >= fit$itmax && stage.at.max == 0)
			stage.at.max <- 1
		if (!is.null(fit[["lastLAObjDiff"]]) &&
			fit$lastLAObjDiff > tol && stage.increased == 0)
			stage.increased <- 1

		LA.at.max <- LA.at.max + stage.at.max
		LA.increased <- LA.increased + stage.increased
	}

	LA.ok <- LA.at.max == 0 && LA.increased == 0
	ok <- optimizer.ok && LA.ok

	cat("Convergence:", if (ok) "OK" else "NOT OK", "\n")
	if (optimized)
		cat(paste0(optimizer.lines, collapse = "\n"), "\n")
	else
		cat("Hyperparameters: not optimized\n")

	if (LA.ok) {
		cat("Laplace approximations: OK\n")
	} else {
		problems <- character()
		if (LA.at.max > 0)
			problems <- c(problems,
				paste(LA.at.max, "reached the maximum number of iterations"))
		if (LA.increased > 0)
			problems <- c(problems,
				paste(LA.increased, "ended with an objective increase"))
		cat("Laplace approximations: NOT OK -",
			paste(problems, collapse = "; "), "\n")
	}

	invisible(ok)
}
