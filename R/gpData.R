

#' gpData constructor
#'
#' Bundle the data for the gaussian process model - construct a gpData object from a named list of tables (data.frame's).
#'
#' For computation efficiency reasons, this package doesn't work in an usual way where a single data.frame which is joined from multiple
#' data.frames. Rather, data is provided as a list of data.frames (we also call them "tables"), where each table can differ in the type of records
#' that each row represents (and thus also in the number of rows).
#' 
#' There are two main cases:
#' 
#' 1. All tables describe the same record type - there is one-to-one correspondence between their rows, and thus they all have the same number of rows.
#' 
#' 2. There are tables representing different record types. In that case, there must be a main table (which is first on the list),
#' and this main table contains so called "grouping factors" (columns). These grouping factors represent all the other record types
#' present in the dataset. 
#' 
#'     - Tables whose rows represent given grouping factor must carry an attribute \code{"fact"} naming that grouping 
#' factor - and they need to have a column of that name, which uniquely identifies its rows. Multiple tables can share the same grouping factor,
#' in which case, there must be one-to-one correspondence between their rows.
#' 
#'     - Tables without grouping factor must represent the same record type as the main table, having one-to-one correspondence between their rows.
#' Throughout the package, absence of a grouping factor is marked as belonging to a special factor \code{"1"}. However, it doesn't need to be denoted with the \code{"fact"} attribute
#' in the input data.
#'
#' Note that if there is only one grouping factor in the data, the main table is not required.
#' That is, all tables can carry the same attribute \code{"fact"} (same grouping factor) without main table to be present.
#' 
#' The function performs validation checks, not tied to any particular formula or object of class `gp`.
#' 
#' @param x Named list of tables (data.frame's). Tables may carry an attribute \code{"fact"} naming a grouping
#' factor column; the column must exist and be a unique identification of the table's rows. See Details below.
#'
#' @return An object of class \code{gpData}
#'
#' @details
#' - Since the special factor value \code{"1"} denotes absence of a grouping factor throughout the package,
#'   naming a factor \code{"1"} is not allowed.
#' - Tables that share the same grouping factor (including the special factor \code{"1"}) must have identical row counts
#'   and there is 1:1 correspondence between their rows.
#'
#' @export
#'

# moje puvodni dokumentace:
#
# This data doesn't have to be standardized, this will be done internally (by the XXX() function).
#
# input:
# - list of tables - data.frame's
#		- if the rows of a given table correspond to some grouping factor, the factor's name
#		  must be specified via attr "fact", and column of the same name must exist in that table - that column
#		  is basically a primary key of that table (in that case, it cannot be a matrix of course)
# - if there tables with different grouping factors (table without grouping factor is considered as a special "grouping factor" for that purpose),
#	the first table is a main table (data.frame); otherwise, no main data.frame needed
# - main df obsahuje vsechny ty grouping factors (ano jsou grouping)
#	- tj. jednoznacne urcuje vsechny ty factors (a tedy typicky bude alespon tolik radek jako ty ostatni df's)
#	- samozrejme ze sam nemuze mit grouping factor
#	- ne nutne urcuje dimenzi GP (dimenzi cov matrix K) - viz case czech atlas without bigm - dimenzi GP urcuje spis ta formula - konkretne ty tables, ktery jsou bez faktoru (ty by mely bejt vsechny stejne velky)
# - all tables with the same grouping factor must have the same number of rows and there must be 1-1 correspondence between the rows. (All tables without grouping factor are considered having the "same grouping factor" for that matter)
#
# !NOTE: dusledky pro dimenzi GP viz komentar v K_matrix.R
#
# CONSEQUENCES OF THE ABOVE CONDITIONS:
# - bud maji vsechny tables radky se stejnym grouping factorem (no grouping factor being special value here), a je 1:1 korespondence mezi radky vsech tabulek
# - nebo se mohou lisit; pak musi byt definovana jedna main table, ktera jednoznacne urcuje vsechny grouping factors (tj. klice vsech ostatnich tabulek)
#	- z toho plyne, ze tato main table musi mit >= radku nez vsechny ostatni tables



# Vnitrni funkcionalita:
# - co je potreba udelat

# samostatne - v constructoru?:
#D - pro kazdy faktor (nebo 1/NULL) zkontrolovat
#D	- ze je pro vsechny tables unique nrow
#D	- a pro ty s fact ze jsou vzdy stejny table[fact,]
#D - z tabulek s fact odstranit fact (jakoby to id)
#D - v main table:
#D 	- zkontrolovat ze jsou vsechny fact
#D	- vytvorit pro kazdy fact fact_idx


# v kontextu modelu - na to bude(ou) dalsi funkce:
# D- zkontrolovat ze data odpovidaji dataReq pro cov
# D- prevest na matrix - ? asi jen to co je potreba pro cov? ANO!
# D	- u tech tabulek co se prevadeji: warnings kdyz int nebo factor (nebo jen obecne kdyz neni numeric?)
# - scaleovat - jen to co je potreba pro cov


# 14.5. vecer: mozna jsem udelal chybu v designu, mozna jsem mel ty data parsovat az podle data requirements z te formule.
# To by melo tu vyhodu, ze by uzivatel nemusel zadavat ten attribut "fact" - to mi doslo, kdyz jsem si predstavil jak pisu manual k gpData
# a musim tam vysvetlovat, co je to grouping factor. Mnohem lepe se to vysvetluje u te formule, v jejim kontextu.
# Uzivatel by tak predal formuli a obycejny list of data.frame's a nemusel by se s tim moc patlat. Ja jsem se rozhodl pro tenhle design
# (delat to samostatne nezavisle na formuli) ze dvou duvodu:
# 1) protoze v nekterych pripadech tam budou tabulky s factorem ktere nebudou ve formula (asi `nests` u cejek)
# 2) protoze si rikam ze pro predict() se bude nejspis hodit vytvaret si predikcni datasety nezavisle na nejakem modelu a jeho formula
# (ale mozna se mylim?)

# todo: zmenit ten check na zacatku a rict is.data.frame() nebo is.matrix a hotovo.

gpData <- function(x)
{
	factors <- list()
	stopifnot(is.list(x))
	first <- TRUE
	for (name in names(x)) {
		d <- x[[name]]
		if (!(is.data.frame(d))) {
			warning("List item `", name, "` - only data.frame's are allowed.")
		}
		fact <- attr(d, "fact")
		if (!is.null(fact)) {
			if (fact == "1")
				stop("Grouping factor cannot be named \"1\".") # because we use it internally in this function for absence of a grouping factor, see fact <- "1" below.
			if (!fact %in% colnames(d))
				stop("List item `", name, "` must contain the grouping factor `", fact, "` as a column.")
			if (!is.data.frame(d))
				stop("Because the list item `", name, "` has a grouping factor, it must be a data.frame.")
			if (any(duplicated(d[[fact]]))) # using this instead of just d[,fact], because with tibble, I would have to do d[,fact,drop=TRUE]...
				stop("List item `", name, "`: the grouping factor `", fact, "` is not unique.")
		}
		if (is.null(fact)) { # no grouping factor here
			fact <- "1"
			id <- 1:nrow(d) # just dummy vector for the comparison below to work as well
		} else {
			id <- d[[fact]]

			# remove the grouping factor (id) from the table (put it in the rownames first)
			rownames(x[[name]]) <- id
			x[[name]][,fact] <- NULL
		}
		# CAREFUL: now fact == 1 instead of NULL

		# test that tables sharing the same factor have the same number of rows and the same fact (id) vector
		if (!fact %in% names(factors)) { # first table with this grouping factor
			factors[[fact]] <- list(
				nrow = nrow(d),
				id = id,
				tables = name
			)
		} else { # another table with this grouping factor - compare with already stored factor - must be the same
			if (factors[[fact]]$nrow != nrow(d))
				stop("List item `", name, "`: tables with the same grouping factor (`", fact, "`) must have the same number of rows.")
			if (!all(factors[[fact]]$id == id))
				stop("List item `", name, "`: tables with the same grouping factor (`", fact, "`) must have the rows in the exact same order (with respect to the id, i.e. the grouping factor)")
			factors[[fact]]$tables <- c(factors[[fact]]$tables, name)
		}
		first <- FALSE
	}
	all.factors <- factors # possibly including "1" for tables with no grouping factor
	factors[["1"]] <- NULL # now, here keep only the "true" grouping factors
	main_name <- NULL
	if (length(all.factors) >= 2) { # if there two tables with different grouping factors
		main <- x[[1]] # this must be the main table
		main_name <- names(x)[1]
		if (!is.data.frame(main))
			stop("The first table in the list, `", main_name, "`, which is supposed to be a main table, should be a data.frame.")
		if (!is.null(attr(main, "fact")))
			stop("The first table in the list, `", main_name, "`, which is supposed to be a main table, cannot have a grouping factor.")
		missing_in_main <- setdiff(names(factors), colnames(main))
		if (length(missing_in_main) > 0) {
			#print(str(missing_in_main))
			#xxx <<- missing_in_main
			stop("The first table in the list, `", main_name, "`, which is supposed to be a main table, is missing these grouping factors: ",
				paste(missing_in_main, collapse = ", "))
		}

		# for each factor, create a corresponding _idx column in the main table
		for (f in names(factors)) {
			idx <- match(main[[f]], factors[[f]]$id) # drop = TRUE would be needed for tibble objects! So I use main[[f]] to make sure it is just a vector
			if (any(is.na(idx)))
				stop("Factor `", f, "` in the main table `", main_name, "` contains values not present in table(s) ", paste(factors[[f]]$tables, collapse = ", "))
			if (length(setdiff(factors[[f]]$id, main[[f]])) > 0)
				warning("Table(s) ", paste(factors[[f]]$tables, collapse = ", "), " contain some values not present in the main table `", main_name, "` (waste of resources).")
			idx_name <- paste0(f, "_idx")
			if (idx_name %in% colnames(main))
				stop("Column named `", idx_name, "` is already present in the main table `", main_name, "` => cannot create index.")
			# create the index now
			x[[main_name]][,idx_name] <- idx
		}
	}
	# double-check the consequences of our validation
	# 1) if there is no main table, all tables must have the same number of rows (we don't check the exact 1-1 correspondence here again, we've already done that)
	if (is.null(main_name))
		stopifnot(all(nrow(x[[1]]) == vapply(x, nrow, integer(1))))
	# 2) main table must have number of rows >= nrows of all other tables
	else
		stopifnot(all(nrow(x[[1]]) >= vapply(x, nrow, integer(1))))
	class(x) <- "gpData"
	attr(x, "factors") <- factors
	attr(x, "main_table") <- main_name # simplify checking for it in further processing
	return(x)
}

#' Subset gpData along given factor and indices
#' 
#' Subsets the gpData object along the given grouping factor `fact` (or along the main table if fact = "1"),
#' taking only the rows/levels with given `indices`.
#' @param x object of class gpData. Does not matter if the data is scaled (gpDataPrepare()'d) or not,
#' 			it works on both.
#' @param fact character. Grouping factor name along which to index, or "1" for no grouping factor (the default), 
#'		which means indexing along the main table (or all tables, if there is no main table; since in the absence of main table 
#' 		all tables must have one-to-one correspondence between their rows).
#' @param ind indices along given factor \code{fact} to subset from the dataset.  If \code{fact = "1"}, indices correspond to the 
#'   main table (or simply all tables, if no main table exists).
#' @return gpData object, subsetted along the given factor and indices.
#' @export
gpDataSubset <- function(x, fact = "1", ind)
{
	stopifnot(max(ind) <= gpDataSize(x, fact))
	stopifnot(is.character(fact))
	if (gpDataHasMainTable(x))
		stopifnot(attr(x, "main_table") == names(x)[1])

	if (fact != "1") {
		stopifnot(!is.null(attr(x, "factors")))
		stopifnot(!is.null(attr(x, "factors")[[fact]]))
		sel_ids <- attr(x, "factors")[[fact]]$id[ind]
		if (gpDataHasMainTable(x))
			ind <- x[[1]][[fact]] %in% sel_ids # get the index along main table
	}
	# subset the main table
	if (gpDataHasMainTable(x)) {
		x[[1]] <- x[[1]][ind,,drop = FALSE]
	}
	factors <- attr(x, "factors")
	# redo all the indices
	f_ind <- list()
	for (fact in names(factors)) {
		fact_idx_name <- paste0(fact, "_idx")
		if (gpDataHasMainTable(x))
			fact_ind <- sort(unique(x[[1]][[fact_idx_name]])) # get the subset index along the factor table
				# done on already subsetted main table; so these are also already subsetted
		else
			fact_ind <- ind
		factors[[fact]]$id <- factors[[fact]]$id[fact_ind]
		factors[[fact]]$nrow <- length(factors[[fact]]$id)
		f_ind[[fact]] <- fact_ind # save it, it will come in handy for the tables!
		if (gpDataHasMainTable(x))
			x[[1]][[fact_idx_name]] <- match(x[[1]][[fact]], factors[[fact]]$id) # reindex the factor in main table
	}
	attr(x, "factors") <- factors
	# redo all the tables, except for the main table
	tables <- names(x)
	if (gpDataHasMainTable(x))
		tables <- names(x)[-1]
	for (tbl in tables) {
		fact <- attr(x[[tbl]], "fact")
		if (is.null(fact))
			fact_ind <- ind # no grouping factor
		else
			fact_ind <- f_ind[[fact]]
		x[[tbl]] <- x[[tbl]][fact_ind,,drop = FALSE]
	}
	return(x)
}





#' internal function
#'
#' check, whether gpData fulfills the requirements given by the gp object (gp$dataReq in particular),
#' i.e. whether all tables required by the formula are present, with correct grouping factors
#'
#
# Maybe this function will change, if we decide to change the semantics so that the gpData() constructor already constructs the gpData according
# to the formula and the dataRequirements implied by that. Then, some of these checks would perhaps be done in that new constructor (but maybe not).
# It either returns TRUE or ends with an error.
#
# !!! mozna zde doplnit parametr `components`, aby tim slo testovat i data pro predikce? To by bylo genialni! A kdyz nebude specif, budou to vsechny komponenty
#		- dodatecne: prijde mi to zbytecny
gpDataCheckReq <- function(gp, gpData)
{
	stopifnot(class(gp) == "gp")
	stopifnot(!is.null(gp$dataReq))
	stopifnot(class(gpData) == "gpData")
	tables_w_grouping_factors <- do.call(c, lapply(attr(gpData, "factors"), function (x) x$tables))
	# Go through all required tables and check if they are there, with correct grouping factors
	for (m in names(gp$dataReq$mats)) { # note that not all the tables in gpData have to be in this set! The main table probably won't be.
		if (!m %in% names(gpData))
			stop("Table `", m, "` is missing in the gpData object.")
		if (gp$dataReq$mats[[m]]$fact != "1") { # the table m should have a grouping factor according to the requirements
			if (!gp$dataReq$mats[[m]]$fact %in% names(attr(gpData, "factors")))
				stop("Factor `", gp$dataReq$mats[[m]]$fact, "` is not present in gpData.")
			if (!m %in% attr(gpData, "factors")[[gp$dataReq$mats[[m]]$fact]]$tables)
				stop("Missing grouping factor `", gp$dataReq$mats[[m]]$fact, "` in table `", m, "`.")
		} else {  # The table m doesn't have a grouping factor
			if (m %in% tables_w_grouping_factors)
				stop("Table `", m, "` shouldn't have any grouping factor, but it does.")
		}
	}
	# Check if all grouping factors are present in the main data.frame
	if (gpDataHasMainTable(gpData) && !all(gp$dataReq$factors %in% colnames(gpData[[1]]))) # tento if funguje i kdyz tam nebudou zadny grouping factors - otestovano!
		stop("These factors are missing in the main data.frame `", names(xx)[1], "`: ",
			paste(setdiff(gp$dataReq$factors, colnames(gpData[[1]])), collapse = ", "))
	return(TRUE)
}

#' Internal function to check whether gpData has a main table
#' @param gpData object of class gpData
#' @return TRUE/FALSE

gpDataHasMainTable <- function (gpData)
{
	stopifnot(class(gpData) == "gpData")
	!is.null(attr(gpData, "main_table"))
}

#' Check whether gpData has been scaled already (gpDataPrepare() has been called on it)
#' 
#' @param gpData object of class gpData
#' @return TRUE/FALSE
#' @export
gpDataIsScaled <- function (gpData)
{
	stopifnot(class(gpData) == "gpData")
	!is.null(attr(gpData, "gpDataPrepared"))
}

# Prepares the data for the model:
#	- converts tables (data.frames) to matrices
#	- scales the matrices (to mean = 0 and sd = 1) that need it (scale = TRUE in cov_funcs.R)
#		- if the gp object already has scaled training data (gp$data), this function will take scaling from there, and
# It issues an error, if there are any non-numeric columns,
#'@export
gpDataPrepare <- function(gp, gpData)
{
	stopifnot(class(gp) == "gp")
	stopifnot(!is.null(gp$covComp))
	stopifnot(class(gpData) == "gpData")
	if (!is.null(gp[["data"]]))
		stopifnot(all(names(gpData) %in% names(gp[["data"]]))) # if training dataset is already present, make sure all tables here were also present in the training
	# note: if training dataset is not present, i.e. gpData is supposed to be the training dataset, we don't have to check if all tables from dataReq are present
	# - this was already done by gpDataCheckReq()
	mats <- intersect(names(gpData), names(gp$dataReq$mats)) # pick all tables that are provided in the data and at the same time used by the formula
	for (m in mats) {
		# convert m to matrix - but first do some checks
		stopifnot(is.data.frame(gpData[[m]]))
		non_num <- !sapply(gpData[[m]], is.numeric)
		if (any(non_num))
			stop("Non-numeric columns in table `", m, "`: ", paste(colnames(gpData[[m]])[non_num], collapse = ", "))
		#int <- sapply(gpData[[m]], is.integer)
		#if (any(int))
		#	warning("Integer columns in table `", m, "`: ", paste(colnames(gpData[[m]])[int], collapse = ", "))
		# convert to matrix now
		if (!gp$dataReq$mats[[m]]$scaling)
			gpData[[m]] <- as.matrix(gpData[[m]])
		else { # scale
			if (is.null(gp[["data"]]))  # training dataset not yet present
										# have to use gp[["data"]] instead of gp$data due to the damn partial matching
				gpData[[m]] <- scale(as.matrix(gpData[[m]]))
			else 					# take scaling from the training dataset
				gpData[[m]] <- scale(as.matrix(gpData[[m]]), center = attr(gp$data[[m]], "scaled:center"), scale = attr(gp$data[[m]], "scaled:scale"))
		}
	}
	attr(gpData, "gpDataPrepared") <- TRUE
	return(gpData)
}


#' Size of the GP dataset (number of rows) along given grouping factor.
#'
#' If the grouping factor \code{fact} is specified, the size unit is this grouping factor, i.e.
#' the function returns number of levels of this grouping factor in the data. 
#' If there is no grouping factor specified (fact = "1", the default), 
#' the dimension is given by the main table, or, if the main table is not present, simply by 
#' all tables (since if the main table is not present, all tables necessarily have the same number 
#' of rows).
#'
#' @param gpData object of class gpData
#' @param fact character. Grouping factor name, or "1" for no grouping factor (the default).
#' @return Integer. Size of the dataset along the given grouping factor.
#' @export
#'
gpDataSize <- function(gpData, fact = "1")
{
	if (fact == "1")
		return(nrow(gpData[[1]])) # which doesn't have to be main table, but if there is no main table, 
								  # all tables have the same number of rows anyway (guaranteed by gpData() constructor)

	# fact is defined
	if (!fact %in% names(attr(gpData, "factors")))
		stop("The fact = `", fact, "` was not found within the grouping factors of this gpData dataset")

	attr(gpData, "factors")[[fact]]$nrow
}




#' Determine the size of the Gaussian Process
#'
#' Determine the size of the Gaussian Process, i.e. the dimension of the covariance matrix, as well as the corresponding grouping factor which defines this size
#' (or "1" if there is no such grouping factor).
#'
#' @details
#' 1. If all covariance components (except intercept) use the same grouping factor (true grouping factor, not "1"), then the dimension of the Gaussian Process is determined
#' by this factor, and the size is given by the number of levels of this factor.
#'
#' 2. In all other cases (no grouping factors, or some components with and some without a grouping factor, or multiple grouping factors used in different components),
#' the size is given by the number of rows of the main table (if the main table is not present, all tables necessarily have the same number of rows). The returned factor is "1",
#' meaning "no grouping factor".
#'
#' @param gp object of class gp
#' @return A list with two elements:
#' \describe{
#' \item{size}{integer - size of the Gaussian Process (dimension of the covariance matrix; i.e. both number of rows and columns, since
#' the covariance matrix is a square matrix)}
#' \item{fact}{character - grouping factor corresponding to the dimension of the Gaussian process which defines the size; or "1" meaning "no grouping factor"}
#' }
#' @export
#'
# dimenze GP (vsechny komponenty zvoleny):
# - pokud tam nejsou grouping factors, zadna table neni main, vsechny pouzite tables have same number of rows a to bude dimenze GP
# - pokud tam jsou grouping factors:
# 		- pokud vsechny komponenty formule (krom interceptu) pouzivaji table s grouping factory, a jedna se o jeden jediny faktor, bude dimenze GP odpovidat tomuto faktoru
#		- ve vsech ostatnich pripadech je to dimenze main table. To se tedy tyka i pripadu kdy zadna komponenta formule neodpovida main table, tj. vsechny krom interceptu maji nejaky
#		  faktor, ale jsou tam alespon dva ruzne faktory - potom je potreba main table na to, aby se tyto propojily
#			!!!! zde by slo zavest hierarchy of the grouping factors a misto toho main by se ty dva (ci vice) faktory mohly propojit v jejich nejblizsim spolecnem parentovi
#				(coz muze byt i jeden z tech dvou faktoru samotnych), ale na to ted pecu

# 2025-10-18: ?? nevrati to presne totez co gpDataSize(gp$data, gp$GP_factor), jak to pouzivam v K_matrix()?
#				asi jo, ale je to cirkularni definice, protoze to gp$GP_factor se inicializuje na zaklade teto funkce!
#		- tj. prejmenovat ji spis na gpDetermineSize() a prezentovat spis jako one time function?
gpDetermineSize <- function (gp)
{
	data <- gp$data
	comps <- gp$covComp_df
	comps <- comps[comps$cov_fun != "1",] # this will drop the intercept component

	factors <- setdiff(na.omit(unique(comps$fact)), "1")
		# get all "true" factors (no NAs, no "1"s)
		# intentionally getting them from the component table rather than from gp$dataReq, just in case I want to generalize it later for particular component selection
	if (length(factors) == 1 && all(!is.na(comps$fact) & comps$fact == factors[1])) {
		# vsechny komponenty formule krom interceptu pouzivaji tento jeden faktor
		size <- attr(data, "factors")[[factors[1]]]$nrow
		fact <- factors[1]
	} else {
		size <- nrow(data[[1]])
		fact <- "1"
	}
	stopifnot(!is.na(fact)) # paranoidni checky, shouldn't be needed, jen na demonstraci ze to splnuje ty predpoklady
	stopifnot(!is.null(fact))
	stopifnot(is.character(fact))
	stopifnot(fact != "")
	list(size = size, fact = fact)
}
