#'@title cc
#'@description Very similar \code{cancor()} of package \code{stats}, but 1) it allows for X and Y to be rank deficient, 2) it allows for categorical variables and 3) it allows for covariates, 4) more (see below).
#'@param X See \code{x} in \code{\link[stats]{cancor}} for a proper documentation.
#'@param Y See \code{y} in \code{\link[stats]{cancor}} for a proper documentation.
#'@param Zx (\code{=NULL} by default) covariates (i.e. nuisance) of \code{X}. If different from \code{NULL}, the \code{X} are residualized by \code{Zx} before entering \code{cc()}. IMPORTANT: if Zx is not NULL, include the intercept (if appropriate!)
#'@param Zy (\code{=Zx} by default) covariates of \code{Y}. Same use of \code{Zx}.
#'@param fill.na replace \code{NA} in \code{X} and \code{Y} with column mean before enter \code{cc()}.
#'@param numb_cc number of (pairs of) canonical correlations to be extracted 
#'@return It returns an \code{acca} object. This object contains the same list as returned by \code{\link[CCA]{cc}}, while it also contains \code{prop_expl_var} which is a \code{list} reporting the proportion of explained (total) variance of \code{X} and \code{Y} by each component (i.e. mode).
#'
#'
#'L=matrix(rnorm(10),10,1)
#'X=matrix(rnorm(50),10,5)
#'Y=matrix(rnorm(70),10,7)
#'Z=matrix(rnorm(20),10,2)
#'
#'X[,1]=X[,1]+2*L
#'Y[,1]=Y[,1]+2*L
#'mod=cc(X,Y,Z)
#'mod
#'
#'ccbiplot(mod)
#'
#'@export cc

cc <- function (X,Y,Zx=NULL,Zy=Zx,numb_cc=NULL,fill.na=FALSE) 
{
  Y=convert2dummies(Y)
  Y=as_named_matrix(Y,"Y")
  Y=fillnas(Y)
  Y <- as.matrix(Y)
  
  X=convert2dummies(X)
  X=as_named_matrix(X,"X")
  X=fillnas(X)
  X <- as.matrix(X)
  
  if(!is.null(Zy))   {
    Zy=convert2dummies(Zy)
    Zy=as_named_matrix(Zy,"Zy")
    Zy=fillnas(Zy)
    Y=residualize(Y,Zy); rm(Zy)
  } else Y=scale(Y,scale=FALSE)
  if(!is.null(Zx))   {
    Zx=convert2dummies(Zx)
    Zx=as_named_matrix(Zx,"Zx")
    Zx=fillnas(Zx)
    X=residualize(X,Zx); rm(Zx)
  }else X=scale(X,scale=FALSE)
  
  Xnames = dimnames(X)[[2]]
  Ynames = dimnames(Y)[[2]]
  ind.names = dimnames(X)[[1]]
  
  
  ###########################
  if ((nr <- nrow(X)) != nrow(Y)) 
    stop("unequal number of rows in 'cancor'")
  ncx <- ncol(X)
  ncy <- ncol(Y)
  if (!nr || !ncx || !ncy) 
    stop("dimension 0 in 'X' or 'Y'")
  ###########################
  qx <- qr(X)
  qy <- qr(Y)
  dx <- qx$rank
  if (!dx) 
    stop("'X' has rank 0")
  dy <- qy$rank
  if (!dy) 
    stop("'Y' has rank 0")
  
  if(dx<ncol(X)){
    svx=svd(X,nu = dx,nv = dx)
    X=svx$u
    svx$u=NULL
    qx <- qr(X)
  } else svx=NULL
  if(dy<ncol(Y)){
    svy=svd(Y,nu = dx,nv = dx)
    Y=svy$u
    svy$u=NULL
    qy <- qr(Y)
  } else svy=NULL
  
  numb_cc=min(numb_cc,dx,dy)
  
  ###################
  mod <- .cc_core(qx,qy,numb_cc)
  ####################
  # compute coeff
  if(TRUE){
    mod$xcoef <- backsolve((qx$qr)[1L:dx, 1L:dx, drop = FALSE], mod$u[,1:numb_cc,drop = FALSE])
    rownames(mod$xcoef) <- colnames(X)[qx$pivot][1L:dx]
    mod$ycoef <- backsolve((qy$qr)[1L:dy, 1L:dy, drop = FALSE], mod$v[,1:numb_cc,drop = FALSE])
    rownames(mod$ycoef) <- colnames(Y)[qy$pivot][1L:dy]
    colnames(mod$xcoef)=paste0("Cx",1:ncol(mod$xcoef))
    colnames(mod$ycoef)=paste0("Cy",1:ncol(mod$ycoef))
  }
  mod$u <- mod$v <- NULL

  if(TRUE) mod$data=list(X=X,Y=Y)
  
  if(TRUE){
    mod=.compute_stats(mod,svx,svy)
  }
  
  mod$call$cc=match.call()
  class(mod) <- c("acca", class(mod))
  return(mod)
}

