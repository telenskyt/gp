
# code to test the function for specific warning suppression

library(Matrix)

rW0 <- Diagonal(n = 11778, x = 1)

inspect_warnings(rW <- as.matrix(rW0), stack = "full")

inspect_warnings(rW <- as.matrix(rW0))

suppress_warnings_from(rW <- as.matrix(rW0), "sparse->dense coercion: allocating vector of size", fun = "asMethod", package = "Matrix")



f <- function()
{
	warning("f fun warns")
}


g <- function()
{
	warning("g fun warns")
	warning("bla bla!")
	f()
}

g()

inspect_warnings(g(), stack = "full")

inspect_warnings(g())

suppress_warnings_from(g(), "fun warns")
suppress_warnings_from(g(), "fun warns", fun = "f")

suppress_warnings_from(g())






