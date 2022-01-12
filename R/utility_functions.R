#' @title Methods for acca objects
#'
#' @description Methods for \code{acca} objects. 
#' The following are methods to extract and manipulate relevant information from
#' a \code{acca} object.
#' 
#' @name acca-method
#' @docType methods

NULL



#' print.acca print method for a acca object.
#' @param x an \code{acca} object
#' @method print acca
#' @docType methods
#' @rdname acca-method
#' @export


print.acca <- function(object, ...) {
  cat("a Canonical Correlation Analysis")
  cat("
              n obs: ", nrow(object$data$X))
  cat("
      n X variables: ", ncol(object$data$X))
  cat("
      n Y variables: ", ncol(object$data$Y))
  cat("
       Correlations: ",sprintf("%.3f",object$cor))
  if(!is.null(object$p_values)){
  cat("
           p-values: ",sprintf("%.3f",object$p_values))
  cat("
        n perms (B): ",attributes(object$p_values)$B)
    
  }
}

#' summary.acca summary method for a acca object.
#' @rdname acca-method
#' @param object an \code{acca}
#' @param ... additional arguments to be passed
#' @method  summary acca
#' @docType methods
#' @export

summary.acca <- function (object, ...) {
  print.acca(object)
}
