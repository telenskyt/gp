

# derivation of cov.SE covariance matrix w.r.t l[i]
# i.e. gradient for each parameter l[i]: dK/dl[i] (this is one iter. of what cov.SE.d1() computes, but optimized memory- and cpu- wise)
# i - number of covariate
# this derivation is:
# d/dl[i] k'(x[j,i], x[k,i]|l[i]) = 1/(2*l[i])*(x[j,i] - x[k,i])^2
# effective computation using one matrix multiplication
#
# !!! special implementation! Takes K for optimization purposes
cov.SE.d1 <- function (K, x1, ls, der_wrt_i)
{
	i <- der_wrt_i	
	K * (cbind(x1[,i], x1[,i]^2, (1 / ls[i] ^ 3)) %*% rbind(-2*x1[,i]/ls[i]^3, (1 / ls[i] ^ 3), x1[,i]^2)) # ls[i]^2 - quadratic here because this is first derivation, otherwise would be ls[i]	
	# orig when ls was already squared: K * (cbind(x1[,i], x1[,i]^2, (1 / ls[i] ^ 2)/2) %*% rbind(-x1[,i]/ls[i]^2, (1 / ls[i] ^ 2)/2, x1[,i]^2)) # ls[i]^2 - quadratic here because this is first derivation, otherwise would be ls[i]
}


#orig function; computes all, too much memory greedy
# not used currently
NOT_USED_cov.SE.d1 <- function (x1, e = NULL, ls) {
  # get gradients (matrices) of the kernel wrt. the parameters
  # CURRENTLY IGNORES e!!
  
  # number of parameters
  n <- length(ls)
  
  # assign vector for gradients
  grads <- list()
  
  # get full covariance matrix
  K <- cov.SE(x1 = x1, e1 = e, ls = ls)
  
  # loop through them
  for (i in 1:n) {
    
    # squared distances
    d2_i <- as.matrix(dist(x1[, i]) ^ 2)
    
    # gradient for each parameter
    grads[[i]] <- K * (1 / ls[i] ^ 2) * d2_i / 2
    
  }
  
  # return as a list
  return (grads)
  
}