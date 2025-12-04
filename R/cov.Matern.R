# cov.Matern(x1, x2, ds) ; v = 3/2
#
# x1 - matrix of covariates for the training dataset (rows - data records, cols - covariates)
# x2 - matrix of covariates for the predicted dataset (rows - data records, cols - covariates)
#
# ds - distance-scale (~ length-scale as defined by Rasmussen & Williams 2006, which are length-scales at the scale of the covariates - 
#		it's not squared length-scale as in cov.SE function!)
# - now we asume it will be a single number (then it is length-scale for all dimensions => length-scale of the distance as well)
# 	- but this function would work with a vector (tength-scales for each dimension), too
#
# argument variants:
# 1) no x2
# 3) with x2

# output: 
#	- covariance matrix of dimensions nrow(x1) x nrow(x1) (variant 1)
#	- covariance matrix of dimensions nrow(x1) x nrow(x2) (variant 3)

#source("my_dist.R")


# at distance ds the covariance is 0.4833577 (at diagonal = 1)
cov.Matern <- function(x1, x2 = NULL, ds) 
{
	if (length(ds) == 1 && ncol(x1) > 1)
		ds <- rep(ds, ncol(x1))
	
	#mstart()
	#dist <- as.matrix(dist(x1 / rep(ds, each = nrow(x1)))) # dist() neni rychlejsi! :)
	dist <- sqrt(my_squared_dist(x1, x2, ds))
	#Dm <<- dist
	#mstop()
    K <- (1 + sqrt(3)*dist) * exp(-sqrt(3)*dist)
}


# first derivation w.r.t ds[i]
# ds - distance-scale (~ length-scale as defined by Rasmussen & Williams 2006, which are length-scales at the scale of the covariates - 
#		it's not squared length-scale as in cov.SE function!)
# - now we asume it will be a single number (then it is length-scale for all dimensions => length-scale of the distance as well)
# i - number of covariate (for the variant with ds vector, which is not supported now)
#
# !! slo by opet tu dist matici uchovavat v nejake K.cache, ale kaslu na to uz. Vypocet teto matice na 4851 x 2 matici x1 trva 1.6 vterin na UZP kompu,
# a navic se bude pocitat jen jedna derivace -> nema to velkou cenu
cov.Matern.d1 <- function (x1, ds, der_wrt = "ds", der_wrt_i = 1)
{
	dist2 <- my_squared_dist(x1 = x1, x2 = NULL, ds)
	dist <- sqrt(dist2)
	# verze kdy ds je vektor:
	# (x <- x1)
	# ddist_dli = 1/(2*dist*ds[i]^2)*(x[j,i] - x[k,i])^2
	#ddist_dli <- (cbind(-x[,i]^2/(2*ds[i]^2), -1/(2*ds[i]^2), x[,i]/ds[i]^2) %*% rbind(1, x[,i]^2, x[,i]))/dist
	# => to plati pro stare ds[i] jako squared distances; kdyz ds[i] jsou ted normalni length-scales, tak je to takto:
	# ddist_dli = 1/(dist*ds[i]^3) * (x[j,i] - x[k,i])^2
	# nebudu to ted predelavat => tuto variantu ted nepodporuji
	#
	# verze kdy ds je jedno cislo pro vsechny dimenze (=> ds pro dist): super simple:
	# ddist_dli <- -dist*/ds
	
	# d cov.Matern / ddist = -3*dist * exp(-sqrt(3)*dist)
	# dK_dli = d cov.Matern / ddist  * ddist_dli ; simplifying expression:
	#dK_dli <- (3/2* exp(-sqrt(3)*dist)/ds[i]^2) * (cbind(x[,i]^2, 1, -2*x[,i]) %*% rbind(1, x[,i]^2, x[,i])) # kdyz ds jsou squared ds
	#dK_dli <- (3*exp(-sqrt(3)*dist)/ds[i]^3) * (cbind(x[,i]^2, 1, -2*x[,i]) %*% rbind(1, x[,i]^2, x[,i])) # kdyz ds jsou normalni ds => tady je asi jeste chyba, ted na to pecu uz
	
	#dK_dli <- -3*dist * exp(-sqrt(3)*dist) * -dist/ds
	if (length(ds) > 1) 
		stop("This variant not implemented at the moment")
	3*dist2/ds * exp(-sqrt(3)*dist)
}
