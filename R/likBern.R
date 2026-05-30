
# bernoulli likelihood; invoked when likelihood is given as string - "bern"
likBern <- function(par, log = FALSE, sum = FALSE)
{	
	getAll(par, warn=FALSE)
	stopifnot(!sum || log) # sum => log
	if (is.matrix(y) && !is.matrix(p) && isTRUE(log) && isTRUE(sum)) {
		# dbinom() would work perfectly anyway for x being matrix and prob vector,
		# but this little optimization - pooling the likelihood across rows - significantly speeds up any RTMB::MakeTape created from this vectorized case! (e.g. in Fisher_sampling)		
		# but it requires that the function calculates the sum already (sum = TRUE)
		y.rs <- rowSums(y)
		stopifnot(length(y.rs) == length(p))
		
		sum(y.rs*log(p) + (ncol(y)-y.rs)*log1p(-p))
			# the only problem here is for p exactly 0 or 1, the result will be 0*log(0) = NaN, see https://chatgpt.com/c/6a19f4be-1e6c-83eb-a785-3c2c2f0bae1d
			# dbinom() has this handled
			# but lets ignore it for now, hopefully it won't happen
	} else {
		res <- dbinom(x = y, size = 1, prob = p, log = log)
		if (isTRUE(log) && isTRUE(sum))
			sum(res)
		else
			res
	}
}


# samples ~ columns
# had to write my own; using RTMB's $simulate() was sketchy (the y <- OBS(y) even converted matrix to vector) and slow
simBern <- function (par, samples)
{
	list(y = matrix(rbinom(n = samples*length(par$p), size = 1, prob = par$p), ncol = samples))
}
