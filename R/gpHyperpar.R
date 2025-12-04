
#source("cov_funcs.R")


# build default hyperparameter table

gpHyperparDefaults <- function (gp)
{


# cc - covariance component
# loop over all cc's

	R <- c()
	for (cc in names(gp$covComp)) {
		ccomp <- gp$covComp[[cc]] # current covariance component
		stopifnot(ccomp$cov_fun %in% names(cov_funcs))
		cov_fun <- cov_funcs[[ccomp$cov_fun]]
		# insert sigma2 of that component (everyone has it)
		R <- rbind(R, tibble(
			component = cc, 
			cov_fun = NA, 
			hyperpar = "sigma2",
			len = 1,
			i = 1,
			table = NA,
			var = NA,
			start = 1, 
			low = 1e-7, # stage 1: 1e-2 !!!
			up = 100,
			value = 1,
			fixed = FALSE,
			optim.link = "log",
			prior = if (ccomp$cov_fun == "1") 
						list(quote(uniform_lp(x))) # for intercept
					else
						list(quote(sigma2_exp_lp(x, lambda = 3)))
			#prior.lambda = 
			#?? related table name?
		))		
		# insert cov_fun-specific hyperparameters
		for (h_name in names(cov_fun$hyperpar)) {
			h <- cov_fun$hyperpar[[h_name]]
			len <- 1 # need to determine manually here
			var <- NA
			stopifnot(!is.null(ccomp$mat) && nchar(ccomp$mat) > 0 && ccomp$mat %in% names(gp$obsdata))
			if (h_name %in% c("ls", "sigma2_int", "sigma2_slope")) {
				# TODO: zautomatizovat len & var, aby se daly definovat v cov_funcs.R
				len <- ncol(gp$obsdata[[ccomp$mat]])
				var <- colnames(gp$obsdata[[ccomp$mat]]) # mozna bude potreba otestovat ze existujou ty colnames? !!!
			} else if (h_name == "sigma2_diag") {
				len <- ncol(gp$obsdata[[ccomp$mat]]) + 1
				var <- c("1", colnames(gp$obsdata[[ccomp$mat]])) # !!! or "(Intercept)"?
					# mozna bude potreba otestovat ze existujou ty colnames? !!!
			}
			R <- rbind(R, tibble(
				component = cc, 
				cov_fun = ccomp$cov_fun, 
				hyperpar = h_name,
				len = len,
				i = 1:len,
				table = ccomp$mat,
				var = var,
				start = h$start, 
				low = h$low,
				up = h$up,
				value = h$start,
				fixed = FALSE,
				optim.link = h$optim.link,
				prior = list(h$prior)
				#prior =
				#prior.lambda = 
				#?? related table name?
			))
		}
	}
	# Make sure the prior is nicely printed (see https://stackoverflow.com/a/79761799/):
	ctl_new_pillar.nice_list_tbl <<- function(controller, x, width, ..., title = NULL) {
		if (!is.list(x)) { # might want to add &!is.data.frame(x) etc.
			return(NextMethod())
		}
		pillar::new_pillar(list(
			#title = pillar::pillar_component(pillar::new_pillar_title(title)), # title - name of the column, just right above the <type>
			title = pillar::new_pillar_title(title), # title - name of the column, just right above the <type>
			type  = pillar::new_pillar_component(list("<list>"), width = 6),
			#data  = pillar::new_pillar_shaft_simple(deparse(x), align = "left") # this is nonsense
			data  = pillar::new_pillar_shaft_simple(x, align = "left") # this works too for class "call"!
			#data  = pillar::new_pillar_shaft_simple(sapply(x, deparse), align = "left") # this works too
		))
	}	
	class(R) <- c("nice_list_tbl", class(R))
	R
}


#!!! starting
#		date.ls = 8,
#		time.ls = 2,
#		!
		
#' @export
# import hyperparameter vector `h`, given by the optimizer, into our hyperparameter table.
gpHyperparImportVector <- function(gp, h, col = "value")
{
	gp$hyperpar[!gp$hyperpar$fixed, col] <- h
	gp$hyperpar
}

#' @export
# export a column from hyperparameter table, for passing it to optimizer (thus skipping the fixed parameters)
gpHyperparExportVector <- function(gp, col = "value")
{
	gp$hyperpar[!gp$hyperpar$fixed, col, drop = TRUE]
}


# get index of i-th optimized (non-fixed) hyperparameter (numbered from 1) within all hyperparameters
# i: index of the hyperparameter within non-fixed (optimized) hyperparameters
#' @export
gpHyperparIdx <- function (gp, i)
{
	stopifnot(i >= 1)
	stopifnot(i <= sum(!gp$hyperpar$fixed))
	which(!gp$hyperpar$fixed)[i]
}

#' @export
# converts a single column of the hyperparameter table to hierarchically structured list of lists of vectors (for each component and hyperparameter)
gpHyperparList <- function (gp, col = "value")
{
	# https://stackoverflow.com/questions/46616791/split-data-frame-by-two-factors
	gp$hyperpar %>% 
		split(~component) %>% 
		map(function (d) split(d[[col]], d$hyperpar))
}

# convert hyperparameter vector h to optimization scale (using link function)
gpLink <- function (gp, h)
{
	res <- c()
	indices <- which(!gp$hyperpar$fixed)
	stopifnot(length(h) == length(indices))
	for (i in 1:length(h)) { # stupid method, but simple
		x <- make.link2(gp$hyperpar[indices[i], "optim.link", drop = TRUE])
		res[i] <- x$linkfun(h[i])
	}
	res
}

# convert optimization scale parameter vector to hyperparameter vector (using inverse link function)
gpLinkInv <- function (gp, par)
{
	res <- c()
	indices <- which(!gp$hyperpar$fixed)
	stopifnot(length(par) == length(indices))
	for (i in 1:length(par)) { # stupid method, but simple
		x <- make.link2(gp$hyperpar[indices[i], "optim.link", drop = TRUE])
		res[i] <- x$linkinv(par[i])
	}
	res
}

# gradient of the inverse link, w.r.t. the optimization scale parameter vector
gpLinkInvDer <- function (gp, par)
{
	res <- c()
	indices <- which(!gp$hyperpar$fixed)
	stopifnot(length(par) == length(indices))
	for (i in 1:length(par)) { # stupid method, but simple
		x <- make.link2(gp$hyperpar[indices[i], "optim.link", drop = TRUE])
		res[i] <- x$linkinvder(par[i])
	}
	res
}

# evaluate prior of the hyperparameter vector h 
# (returns vector)
gpPrior <- function (gp, h)
{
	res <- c()
	indices <- which(!gp$hyperpar$fixed)
	stopifnot(length(h) == length(indices))
	for (i in 1:length(h)) { # stupid method, but simple
		res[i] <- eval(gp$hyperpar$prior[[indices[i]]], list(x = h[i]))
	}
	res	
}

gpPriorGradient <- function (gp, h)
{
	res <- c()
	indices <- which(!gp$hyperpar$fixed)
	stopifnot(length(h) == length(indices))
	for (i in 1:length(h)) { # stupid method, but simple
		call <- gp$hyperpar$prior[[indices[i]]]
		call[[1]] <- sym(paste0(as.character(call[[1]]), "g")) # call_modify() doesn't work for me, so have to do it directly
		res[i] <- eval(call, list(x = h[i]))
	}
	res	
}

#' @export
gpHyperparCheck <- function(gp, h, tol = sqrt(.Machine$double.eps))
{
	stopifnot(length(h) == nrow(gp$hyperpar))
	stopifnot(all(h >= gp$hyperpar$low - tol))
	#stopifnot(all(h >= gp$hyperpar$low))
		#stop("all(h >= gp$hyperpar$low) is not TRUE, see indices ", which(!(h >= gp$hyperpar$low)))
	#stopifnot(all(h <= gp$hyperpar$up))
	stopifnot(all(h <= gp$hyperpar$up + tol))
}

#' Check consistency of the hyperp
#'
#' @export
gpHyperparCheckAll <- function (gp)
{
	gpHyperparCheck(gp, gp$hyperpar$start)
	gpHyperparCheck(gp, gp$hyperpar$value)
}

# get the starting values from another model (gp0) wherever possible
# will set the hyperpar columns `start` and `value` based on `value` from model gp0 wherever it will match the component, hyperparameter name and variable name
#' @export
gpHyperparStartFromModel <- function(gp, gp0)
{
	hy <- gp$hyperpar
	prev_hy <- gp0$hyperpar	
	hy2 <- hy %>% 
		left_join(select(prev_hy, c(component, hyperpar, var, prev_optimum = value)), join_by(component, hyperpar, var)) %>%
		mutate(start = coalesce(prev_optimum, start), value = start, .keep = "unused") 
	gp$hyperpar <- hy2
	gpHyperparCheckAll(gp)
	gp
}



#hyperpar 
#	D- decode/encode 
#	- link functions
#	- priors
#	- range - not needed
#	- starting values  - not needed
	