#' A Canonical Correlation Analysis with Inferential Guaranties
#'
#' @description It performs Canonical Correlation Analysis and provides inferential guaranties on the correlation components. The p-values are computed following the resampling method developed in Winkler, A. M., Renaud, O., Smith, S. M., & Nichols, T. E. (2020). Permutation inference for canonical correlation analysis. NeuroImage, <doi:10.1016/j.neuroimage.2020.117065>. Furthermore, it provides plotting tools to visualize the results.
#' @examples
#' set.seed(1)
#' X=matrix(rnorm(500),100,5)
#' Y=matrix(rnorm(700),100,7)
#' Z=matrix(rnorm(200),100,2)
#' mod=cc(X,Y,Z)
#' mod
#' 
#' ccbiplot(mod)
#' 
#' mod=cc_inference(mod, B = 100, numb_cc = 3)
#' mod
#' 
#' @docType package
#'
#' @author Livio Finos
#' @name acca-package
#' @importFrom methods is
#' @importFrom stats cor model.matrix qchisq rbinom var

NULL