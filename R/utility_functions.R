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
#' @param x a acca object
#' @method print acca
#' @docType methods
#' @rdname acca-method
#' @export


print.acca <- function(x, ...) {
  cat("Canonical Correlation Analysis:")
  cat("
      n observations: ", nrow(x$data$X))
  cat("
      n X variables: ", ncol(x$data$X))
  cat("
      n Y variables: ", ncol(x$data$Y))
  cat("
      Correlations: ",x$cor)
}

#' summary.acca summary method for a acca object.
#' @rdname acca-method
#' @param object a acca object
#' @param ... additional arguments to be passed
#' @method  summary acca
#' @docType methods
#' @export

summary.acca <- function (object, ...) {
  print.acca(object)
}
