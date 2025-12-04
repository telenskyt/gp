 
# wrapper for stats::make.link()
make.link2 <- function(link)
{
	if (link == "cubicrootlog") {
		lnk <- list()
		lnk$name <- link
		lnk$linkinv <- function(x, a = 0.5, c = 0.5) exp(a*x^3+c*x) # from optimization scale to hyperparameter scale
		lnk$linkinvder <- function(x, a = 0.5, c = 0.5) exp(a*x^3+c*x)*(3*a*x^2+c)
		lnk$linkfun <- function (y, a = 0.5, c = 0.5) # from hyperparameter scale to optimization scale
		# when the hyperparameter is a length-scale, I use a = 0.5, c = 0.5
		# when it was a squared length-scale, I used a = 1, c = 1
		{
			# inv. fce te theta.link.inv()
			# roots podle Cardano (kdyz b = 0)
			#y <- 2
			p <- c/a
			q <- -log(y)/a
			cubic.root(-q/2+sqrt(q^2/4+p^3/27))+cubic.root(-q/2-sqrt(q^2/4+p^3/27))
		}
		lnk$mu.eta <- lnk$linkinvder
		lnk$valideta <- function (eta) TRUE # eta ~ x
		class(lnk) <- "link-glm"
	} else {
		lnk <- make.link(link)
		lnk$linkinvder <- lnk$mu.eta # this name really makes much more sense - it is a derivative of inverse link function
	}
	lnk
}
