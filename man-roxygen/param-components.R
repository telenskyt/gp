#' @param components of the covariance formula. One of the following:
#'	- NULL : all components will be selected (the default)
#'	- character vector: components to be selected
#'	- formula: additive formula specifying components, same syntax as in \link[stats:update.formula]{update.formula}, e.g. 
#' \tabular{ll}{
#'   \code{~.}          \verb{   } \tab all components \cr
#'   \code{~.-env-spat} \verb{   } \tab all components except \code{env} and \code{spat} \cr
#'   \code{~spat+year}  \verb{   } \tab only components \code{spat} and \code{year} \cr
#' }
#'
