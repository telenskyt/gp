
# h - vector of hyperparameters that are being optimized (not fixed), optim.link funkci musi resit volajici, zde se uz okolo tohoto nic neresi
# NEBO LIST OF HYPERPARAMETERS??? jako by mel jit do d0, d1, ... i kdyz tam taky zatim nereseno
mnfun <- function (gp, h)
{
	rep(0, gp$GP_size)
}
