#' a canonical correlation analysis 
#'
#' @description a canonical correlation analysis
#' @examples
#' set.seed(1)
#' X=matrix(rnorm(500),100,5)
#' Y=matrix(rnorm(700),100,7)
#' Z=matrix(rnorm(200),100,2)
#' mod=cc(X,Y,Z)
#' mod
#' 
#' ggbiplot2(mod)
#' 
#' mod=cc_inference(mod,nperms = 100,numb_cc = 3)
#' mod
#' 
#' @docType package
#'
#' @author Livio Finos
#' @name acca-package
NULL