
# basic calculus functions, especially with Matrices


#' Tests if matrix M is numeric, with all finite, non-NA values. Works on dense as well as sparse matrices.
#' @export
# source https://chatgpt.com/c/6a157ec4-2d4c-83eb-a76b-f2459bfaf8ef
is_numeric_finite <- function(M) 
{
	if (is(M, "sparseMatrix")) {
		is(M, "dMatrix") && !anyNA(M@x) && !any(is.infinite(M@x))
	} else if (is(M, "Matrix")) {
		is(M, "dMatrix") && !anyNA(M) && !any(is.infinite(M))
	} else {
		is.numeric(M) && !anyNA(M) && !any(is.infinite(M))
	}
}

#' Returns number of structural non-zeros in a sparse matrix
#' @export
# source https://chatgpt.com/c/6a15846d-df40-83eb-8937-1e097e102675
snnzero <- function(x) 
{
	if (!methods::is(x, "sparseMatrix"))
		stop("`x` must inherit from Matrix::sparseMatrix.", call. = FALSE)

	if (methods::is(x, "CsparseMatrix") || methods::is(x, "RsparseMatrix")) {
		p <- methods::slot(x, "p")
		return(unname(p[length(p)]))
	}

	if (methods::is(x, "TsparseMatrix")) {
		return(length(methods::slot(x, "i")))
	}

	if (methods::is(x, "diagonalMatrix")) {
		return(dim(x)[1L])
	}

	if (methods::is(x, "indMatrix")) {
		return(length(methods::slot(x, "perm")))
	}

	stop("Unsupported sparse Matrix class: ", paste(class(x), collapse = ", "))
}

#' Mask dense matrix D to the structural non-zeros of sparse matrix S
#' @param D dense matrix 
#' @param S sparse matrix
#' @returns a sparse matrix, which is a matrix D but masked to the structural non-zeros of matrix S
#' @export
# source https://chatgpt.com/c/6a15aa9b-1bec-83eb-a08c-e17da9badfe5 ; manually improved the suggestion
mask <- function(D, S) 
{
	if (!inherits(S, "sparseMatrix")) {
		stop("S must inherit from 'sparseMatrix'.", call. = FALSE)
	}

	stopifnot(is.matrix(D))
	
	if (!identical(dim(D), dim(S))) {
		stop("dim(D) must equal dim(S).", call. = FALSE)
	}

	P <- as(S, "CsparseMatrix")

	nr <- nrow(D)
	colNum <- rep.int(seq_len(ncol(D)) - 1L, diff(P@p)) # column number, numbered from 0

	methods::new(
		"dgCMatrix",
		Dim      = P@Dim,
		Dimnames = dimnames(D),
		i        = P@i,
		p        = P@p,
		x        = as.numeric(D[P@i + 1L + colNum * nr])
	)
}





mask_too_slow_even_with_MatrixExtra_package_loaded <- function(D, S) 
{
	stopifnot(inherits(S, "sparseMatrix"))

	D <- as.matrix(D)
	if (!identical(dim(D), dim(S))) {
		stop("dim(D) must equal dim(S).", call. = FALSE)
	}

	as(as(S, "nMatrix"), "CsparseMatrix") * D
}

#' Flattens sparse matrix in column-by-column order
#' This flattening procedure must be exactly the same as performed by RTMB Tape class, method $jacfun(sparse = TRUE),
#' see https://github.com/kaskr/RTMB/issues/85
#'
#' @export
flatten <- function(S)
{
	stopifnot(is(S, "CsparseMatrix"))
	S@x
}
