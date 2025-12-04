# moje distancni fce, pomocna pro covariancni matice
# vykuchana z moji verze cov.SE, da se totiz pouzit i pro jine cov. matice i pro jine ucely!
# vrati matici druhych mocnin euklid. vzdaleností, ktere jsou v jednotlivych dimenzich skalovany length-scalemi
# my_squared_dist() je rychlejsi nez dist()! , a umoznuje variantu se dvema maticemi x1, x2
# 	a sqrt(my_squared_dist()) je rychlejsi nez as.matrix(dist()) - ale ne nez dist() sám
#
# my_squared_dist(x1, x2, l)
# x1 - matrix of covariates for the training dataset (rows - data records, cols - covariates)
# x2 - matrix of covariates for the predicted dataset (rows - data records, cols - covariates)
# ls - length-scales = as defined by Rasmussen & Williams 2006, which are length-scales at the scale of the covariates
# original memory requirements: n1*n2*n3
# new memory requirements: n1*n2 (since v1)
#
# changelog:
# v3 - single matrix multiplication with added rows and cols (the trick I applied in gpFitLaplace.R as well); further speedups and memory savings for 1,3,4,5, see cov.SE.test_LOG.txt
#
# argument variants:
# 1) no x2
# 3) with x2
#
# output: 
#	- matrix of dimensions nrow(x1) x nrow(x1) (variant 1)
#	- matrix of dimensions nrow(x1) x nrow(x2) (variant 3)
my_squared_dist <- function(x1, x2 = NULL, ls) 
{
  n1 <- nrow(x1)
  n2 <- ifelse(is.null(x2), n1, nrow(x2))
  n3 <- ncol(x1)
  
#  print(dim(x1))
 # print(n1)
  #print(n2)
  #print(n3)
  
  stopifnot(all(ls > 0)) # not sure about this... test it!

	# trick to big speedup: apply length-scales to the covariates right away!
	x1 <- x1 / rep(ls, each = n1)
	if (is.null(x2)) { # variant 1)
		# if no second matrix do with distance matrices for speed up
		#sumdiffs <- dist(x1)^2  # matrix n1 x n1 (=n2): SQUARED distances, scaled with lengths
			# old time: 5.670s, new time: 2.640s, old mem: 2173.9Mb, new mem: 907.3Mb. Test 1) no x2, no e1/e2, for n3 = 12, n1 = 4873, n2 = 4900
		# this is even faster than the dist() function !!!! Probably because class 'dist' must be transformed to matrix and that is time costly
		x1_2 <- apply(x1^2, 1, sum)
		#sumdiffs <- matrix(x1_2, nrow = n1, ncol = n2) + matrix(x1_2, nrow = n1, ncol = n2, byrow = TRUE) - 2 * x1 %*% t(x1)
			# old time: 5.330s, new time: 1.720s, old mem: 2173.8Mb, new mem: 727.5Mb. Test 1) no x2, no e1/e2, for n3 = 12, n1 = 4873, n2 = 4900
		
		# computing the distances by a single matrix multiplication, as x[i]^2 + x[j]^2 - 2*x[i]*x[j], by putting two extra columns/rows to add the quadratic terms to the sum
		sumdiffs <- cbind(x1, x1_2, 1) %*% rbind(t(-2*x1), 1, x1_2)
			# old time: 4.610s, new time: 1.340s, old mem: 2083.3Mb, new mem: 366.9Mb. Test 1) no x2, no e1/e2, for n3 = 12, n1 = 4873, n2 = 4900
	} else { # variant 3)
	  x2 <- x2 / rep(ls, each = n2)
      #sumdiffs <- 0
      #for (i in 1:n3) { # perhaps this could be vectorized as well, by replacing the %*% with rep() (or by matrix(), see https://stackoverflow.com/a/14927622/684229 !!!)
	  	#dists_i <-   x1[, i] ^ 2 %*% t(rep(1, n2)) +
		#	rep(1, n1) %*% t(x2[, i] ^ 2) - 2 * x1[, i] %*% t(x2[, i]) # matrix n1 x n2: (x1[row, i] - x2[col, i])^2, scaled with lengths
			
		# w/outer product
	  	#dists_i <-   x1[, i] ^ 2 %*% t(rep(1, n2)) +
		#	rep(1, n1) %*% t(x2[, i] ^ 2) - 2 * x1[, i] %o% x2[, i] # matrix n1 x n2: (x1[row, i] - x2[col, i])^2, scaled with lengths	
				# slightly smaller consumption of memory but slightly worse time
		
		#dists_i <- (matrix(x1[,i], nrow = n1, ncol = n2) - matrix(x2[,i], nrow = n1, ncol = n2, byrow = TRUE))^2 # surprisingly slightly longer time and slightly more memory
        #sumdiffs <- sumdiffs + dists_i
      #}
	  # Brutal speed-up! :-) old time: 7.830s, new time: 1.780s, old mem: 3279.1Mb, new mem: 735.0Mb. Test 3) with x2, no e1/e2, for n3 = 12, n1 = 4873, n2 = 4900
	  #sumdiffs <- matrix(apply(x1^2, 1, sum), nrow = n1, ncol = n2) + matrix(apply(x2^2, 1, sum), nrow = n1, ncol = n2, byrow = TRUE) - 2 * x1 %*% t(x2)
	  sumdiffs <- cbind(x1, apply(x1^2, 1, sum), 1) %*% rbind(t(-2*x2), 1, apply(x2^2, 1, sum))
	}	
	sumdiffs[sumdiffs < 0] <- 0 
		# numerical fix, because of Intel oneAPI - oneMKL there are sometimes values like -7.275958e-12!
		# I am not just doing it on the diagonal, because it can happen also anywhere if there is the same location duplicated.
	return(sumdiffs)
}
