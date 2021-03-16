#'cc
#'@description Very similar \code{cc()} of package \code{CCA}, but 1) it allows for X and Y to be rank deficient, 2) it allows for categorical variables and 3) it allows for covariates, 4) more (see below).
#'@param X See \code{\link[CCA]{cc}} for a proper documentation.
#'@param Y See \code{\link[CCA]{cc}} for a proper documentation.
#'@param Zx (\code{=NULL} by default) covariates (i.e. nuisance) of \code{X}. If different from \code{NULL}, the \code{X} are residualized by \code{Zx} before entering \code{cc()}. IMPORTANT: if Zx is not NULL, include the intercept (if appropriate!)
#'@param Zy (\code{=Zx} by default) covariates of \code{Y}. Same use of \code{Zx}.
#'@param fill.na replace \code{NA} in \code{X} and \code{Y} with column mean before enter \code{cc()}.
#'@return The returned list is the same list as returned by \code{\link[CCA]{cc}}, while it also contains \code{prop_expl_var} which is a \code{list} reporting the proportion of explained (total) variance of \code{X} and \code{Y} by each component (i.e. mode).
#'@export
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
#'ggbiplot2(mod)


cc<-function (X,Y,Zx=1,Zy=Zx,numb_cc=NULL,fill.na=FALSE) 
{
  Y=convert2dummies(Y)
  Y=as_named_matrix(Y,"Y")
  Y=fillnas(Y)
  X=convert2dummies(X)
  X=as_named_matrix(X,"X")
  X=fillnas(X)
  
  if(!is.null(Zy))   {
    Zy=convert2dummies(Zy)
    Zy=as_named_matrix(Zy,"Zy")
    Zy=fillnas(Zy)
    Y=residualize(Y,Zy); #rm(Zy)
  } else Y=scale(Y,scale=FALSE)
  if(!is.null(Zx))   {
    Zx=convert2dummies(Zx)
    Zx=as_named_matrix(Zx,"Zx")
    Zx=fillnas(Zx)
    X=residualize(X,Zx); #rm(Zx)
  }else X=scale(X,scale=FALSE)
  
  
  
  Xnames = dimnames(X)[[2]]
  Ynames = dimnames(Y)[[2]]
  ind.names = dimnames(X)[[1]]
  
  
  ###########################
  svx=.svd(X);
  svy=.svd(Y);

  mod=.cc_core(svx,svy,numb_cc = numb_cc)

  rownames(mod$xcoef)=colnames(X)
  rownames(mod$ycoef)=colnames(Y)
  colnames(mod$xcoef)=paste0("Cx",1:ncol(mod$xcoef))
  colnames(mod$ycoef)=paste0("Cy",1:ncol(mod$ycoef))
  mod$data=list(X=X,Y=Y,Zx=Zx,Zy=Zy)
  
  mod=.compute_stats(mod)
  mod$call$cc=match.call()
  class(mod) <- c("acca", class(mod))
  return(mod)
}

