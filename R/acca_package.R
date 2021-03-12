#' a canonical correlation analysis 
#'
#' @description a canonical correlation analysis
#' @examples
#' set.seed(1)
#' X=matrix(rnorm(500),100,5)
#' Y=matrix(rnorm(700),100,7)
#' Z=matrix(rnorm(200),100,2)
#' mod=cc(X,Y,Z)
#' mod$cor
#' mod$scores$corr.X.xscores
#' mod$prop_expl_var
#' ggbiplot2(mod)
#' 
#' @docType package
#'
#' @author Livio Finos
#' @name acca-package
NULL