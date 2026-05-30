
# can be given as string - just "bern"
likBern <- function(par, log = FALSE)
{	
	getAll(par, warn=FALSE)
	dbinom(x = y, size = 1, prob = p, log = log)
}


# samples ~ columns
simBern <- function (par, samples)
{
	list(y = matrix(rbinom(n = samples*length(par$p), size = 1, prob = par$p), ncol = samples))
}
