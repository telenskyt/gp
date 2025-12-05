

#' Validates the components argument - the selection of components
#'
#' 
#' @param gp object of class gp
#' @template param-components
#' 
#' @returns character vector of component names
#' @export 
#' 


#
#		`~.`    		all components
#		`~.-env-spat`	all components except env and spat
#		`~spat+year`    only components spat and year

validate_components <- function (gp, components)
{
	if (is.null(components)) {
		components <- names(gp$covComp) # select all components
		return(components)
	}
	if (is.formula(components)) {
		f0 <- as.formula(paste0("~", paste(names(gp$covComp), collapse = "+"))) # make formula from all components
		f <- update.formula(f0, components) # update it according to a given formula
		if (f[[1]] != "~")
			stop("formula must contain ~")

		if (length(f) != 2) 
			stop("left side of the formula should be empty")	

		f <- f[[2]]
		components <- c() # the converted vector of components
		
		repeat { # over summed terms

			if (length(f) > 1 && f[[1]] == '+') {
				stopifnot(length(f) == 3)
				t <- f[[3]]
				f <- f[[2]]
			} else {
				t <- f
				f <- NULL
			}
			if (!is.symbol(t) || length(t) != 1) {
				stop("The term '", deparse(t), "' should be only a simple component name.")
			}
			components <- c(deparse(t), components) # add before the existing list, so that order is preserved

			if (is.null(f))
				break
		}
		# the validation of the components will continue below		
	}
	stopifnot(is.character(components))
	stopifnot(length(components) > 0)
	invalid <- setdiff(components, names(gp$covComp))
	if (length(invalid) > 0)
		stop("invalid components: ", invalid)
	return(components)
}

# components: see validate_components()
#
# vrati tabulku (df) s prehledem komponent a jejich sizes v training dataset
# !!! uvazoval jsem, ze by bylo asi snazsi, aby tuto tabulku (krom tech sizes) rovnou generovala funkce gpFormula()  - ale to az pripadne na nejaky redesign
gpComponentsTable <- function(gp, data = gp$obsdata, components = NULL)
{
	components <- validate_components(gp, components)
	
	covComp <- gp$covComp
	
	# Helper function to replace NULLs with NAs in a list
	null_to_na <- function(x) {
		lapply(x, function(e) if (is.null(e)) NA else e)
	}

	# Convert covComp to a data frame
	covComp_df <- do.call(rbind, lapply(components, function(name) {
		data.frame(component = name, as.data.frame(null_to_na(covComp[[name]]), stringsAsFactors = FALSE))
	}))
		# 2025-10: null_to_na uz by asi nemelo byt potreba ted? Ale nechavam jeste...

	covComp_df[!is.na(covComp_df$mat), 'nrow'] <- vapply(covComp_df[!is.na(covComp_df$mat), 'mat'] , function (x) nrow(data[[x]]), integer(1))


	#covComp_df[is.na(covComp_df$mat) & !is.na(covComp_df$fact) & covComp_df$fact != "1", 'nrow'] <- vapply(covComp_df[is.na(covComp_df$mat) & !is.na(covComp_df$fact) & covComp_df$fact != "1", 'fact'], function (x) attr(data, "factors")[[x]]$nrow, integer(1)) # should be faster than nrow(), ale pro faktory ktery nejsou id v zadne table to nejde - treba (1|year) 
		
	covComp_df[is.na(covComp_df$mat) & !is.na(covComp_df$fact) & covComp_df$fact != "1", 'nrow'] <- vapply(covComp_df[is.na(covComp_df$mat) & !is.na(covComp_df$fact) & covComp_df$fact != "1", 'fact'], function (x) length(unique(data[[1]][[x]])), integer(1))		

	covComp_df
}

