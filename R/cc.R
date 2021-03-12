#'cc
#'@description Very similar \code{cc()} of package \code{CCA}, but 1) it allows for X and Y to be rank deficient, 2) it allows for categorical variables and 3) it allows for covariates, 4) more (see below).
#'@param X See \code{\link[CCA]{cc}} for a proper documentation.
#'@param Y See \code{\link[CCA]{cc}} for a proper documentation.
#'@param Zx (\code{=NULL} by default) covariates (i.e. nuisance) of \code{X}. If different from \code{NULL}, the \code{X} are residualized by \code{Zx} before entering \code{cc()}.
#'@param Zy (\code{=Zx} by default) covariates of \code{Y}. Same use of \code{Zx}.
#'@param fill.na replace \code{NA} in \code{X} and \code{Y} with column mean before enter \code{cc()}.
#'@return The returned list is the same list as returned by \code{\link[CCA]{cc}}, while it also contains \code{prop_expl_var} which is a \code{list} reporting the proportion of explained (total) variance of \code{X} and \code{Y} by each component (i.e. mode).
#'@export

cc<-function (X,Y,Zx=NULL,Zy=Zx,numb_cc=NULL,fill.na=FALSE) 
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
  }
  if(!is.null(Zx))   {
    Zx=convert2dummies(Zx)
    Zx=as_named_matrix(Zx,"Zx")
    Zx=fillnas(Zx)
    X=residualize(X,Zx); #rm(Zx)
  }
  
  
  
  Xnames = dimnames(X)[[2]]
  Ynames = dimnames(Y)[[2]]
  ind.names = dimnames(X)[[1]]
  
  
  ###########################
  X=scale(X,scale=FALSE)
  svx=svd(X);
  Y=scale(Y,scale=FALSE)
  svy=svd(Y);

  if(any(c(svx$d < 1e-12), (svy$d < 1e-12))){
    np=svx$d>1E-12
    rotX=svx$v[,np]
    svx$v=diag(sum(np))
    svx$u=svx$u[,np]
    svx$d=svx$d[np]
    
    np=svy$d>1E-12
    rotY=svy$v[,np]
    svy$v=diag(sum(np))
    svy$u=svy$u[,np]
    svy$d=svy$d[np]
    
    ###############
    mod=.cc_core(svx,svy,numb_cc = numb_cc)
    mod$names$Xnames=Xnames
    mod$xcoef=rotX%*%mod$xcoef
    
    mod$names$Ynames=Ynames
    mod$ycoef=rotY%*%mod$ycoef
    
    mod$scores$corr.Y.xscores=cor(Y,mod$scores$xscores)
    mod$scores$corr.Y.yscores=cor(Y,mod$scores$yscores)
    mod$scores$corr.X.xscores=cor(X,mod$scores$xscores)
    mod$scores$corr.X.yscores=cor(X,mod$scores$yscores)
  } else {
    mod=.cc_core(svx,svy,numb_cc = numb_cc)
  } 
  
  rownames(mod$xcoef)=colnames(X)
  rownames(mod$ycoef)=colnames(Y)
  colnames(mod$xcoef)=paste0("Cx",1:ncol(mod$xcoef))
  colnames(mod$ycoef)=paste0("Cy",1:ncol(mod$ycoef))
  mod$data=list(X=X,Y=Y,Zx=Zx,Zy=Zy)
  
  mod=.compute_stats(mod)
  class(mod) <- c("acca", class(mod))
  return(mod)
}

