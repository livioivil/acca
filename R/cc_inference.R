#'cc_inference
#'@description For each pairs of components, it computes p-values to test the null hypothesis of no correlation between components. The p-values are computed following the resampling method developed in Winkler, A. M., Renaud, O., Smith, S. M., & Nichols, T. E. (2020). Permutation inference for canonical correlation analysis. NeuroImage, 220, 117065. https://doi.org/10.1016/j.neuroimage.2020.117065.
#'@param mod an \code{acca} object.
#'@param B (\code{100} by default) number of random sampling
#'@param alpha_max stop if p-value > alpha_max (\code{alpha_max=0.5} by default).
#'@param numb_cc stop after computing p-values for the first \code{numb_cc} are computed(\code{numb_cc=NULL} - the default - means compute all).
#'@param resamp_type \code{"sign-flip"} (by default) or \code{"permutation"}
#'@param light If \code{TRUE} the orthogonalization of the residuals of the projection on Z is not performed. For large sample size the two methods tend to overlap. 
#'@return It returns an   \code{acca} object (see \code{\link[acca]{cc}}) with p-values for each pair of the \code{numb_cc} components.
#'
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
#'@export
#'
cc_inference <-  function(mod,B=100, alpha_max=.5,numb_cc=NULL,resamp_type="sign-flip",light=FALSE){
  mod$call$cc_inference=match.call()
  n=nrow(mod$data$X)
  resamp_type=match.arg(resamp_type,c("sign-flip","permutation"))
  if(is.null(numb_cc)) numb_cc=length(mod$cor)
  
  
  ####################
  if(!light) {
    mod=.cc_inference_orthogonal(mod,B, alpha_max,numb_cc,resamp_type)
  } else   if(light) {
    mod=.cc_inference_residuals(mod,B, alpha_max,numb_cc,resamp_type)
  }
  ##################

  return(mod)
}


#############

.cc_inference_orthogonal <- function(mod,B, alpha_max,numb_cc,resamp_type){
  resid_matrix <- function(Zy)
    residualizing_matrix(Zy)$Q
  if(resamp_type=="sign-flip"){
    .permute <- function(X,Qx,nred){
      t(Qx*(1-2*rbinom(nred,1,.5))) %*%X # X=Qx%*%mod$data$X
    }
  } else if(resamp_type=="permutation"){
    .permute <- function(X,Qx,nred){
      t(Qx[sample(nred),]) %*%X # X=Qx%*%mod$data$X
    }
  }
  perm_and_cc=function(X,Y,Qx,Qy,nredx,nredy){
    ccp=.cc_core(.permute(X,Qx,nredx), .permute(Y,Qy,nredy),numb_cc = 0)
    ccp$cor[1]
  }
  
  if(is.null(mod$data$Zx)) {
    Qx=diag(nrow(mod$data$X))
    X=mod$data$X
  } else{
    Qx=resid_matrix(mod$data$Zx)
    X=Qx%*%mod$data$X
  }
  
  if(is.null(mod$data$Zy)) {
    Qy=diag(nrow(mod$data$Y))
    Y=mod$data$Y
  } else{
    Qy=resid_matrix(mod$data$Zy)
    Y=Qy%*%mod$data$Y
  }
  
  nredx=nrow(Qx)
  nredy=nrow(Qy)
  
  
  mod$p_values=rep(1,length(mod$cor))
  attr(mod$p_values,which = "B")=B
  
  
  n_nuisx=min(0,ncol(mod$data$Zx))
  n_nuisy=min(0,ncol(mod$data$Zy))
  Zx=cbind(mod$data$Zx,mod$scores$xscores[,1:numb_cc])
  Zy=cbind(mod$data$Zy,mod$scores$yscores[,1:numb_cc])
  
  for(i in 1:numb_cc){
    mod$p_values[i]=(sum(replicate(B,perm_and_cc(X,Y,Qx,Qy,nredx,nredy))>=mod$cor[i])+1)/(B+1)
    if(mod$p_values[i]>= alpha_max) return(mod)
    Qx=resid_matrix(Zx[,1:n_nuisx+i])
    Qy=resid_matrix(Zy[,1:n_nuisy+i])
    X=Qx%*%mod$data$X
    Y=Qy%*%mod$data$Y
  }
 
  return(mod)
}

#############

.cc_inference_residuals <- function(mod,B, alpha_max,numb_cc,resamp_type){
  n=nrow(mod$data$X)
  if(resamp_type=="sign-flip"){
    .permute <- function(X,n){
      (1-2*rbinom(n,1,.5))*X # X=IH%*%mod$data$X
    }
  } else if(resamp_type=="permutation"){
    .permute <- function(X,n){
      X[sample(n),] # X=IH%*%mod$data$X
    }
  }
  
  perm_and_cc=function(X,Y,n){
    ccp=.cc_core(.permute(X,n), Y,numb_cc =  0)
    ccp$cor[1]
  }
  
  if(is.null(mod$data$Zx)) {
    X=mod$data$X
  } else{
    X=residualize(mod$data$X,mod$data$Zx)
  }
  
  if(is.null(mod$data$Zy)) {
    Y=mod$data$Y
  } else{
    Y=residualize(mod$data$Y,mod$data$Zy)
  }
  
  
  mod$p_values=rep(1,length(mod$cor))
  attr(mod$p_values,which = "B")=B
  
  
  n_nuisx=min(0,ncol(mod$data$Zx))
  n_nuisy=min(0,ncol(mod$data$Zy))
  Zx=cbind(mod$data$Zx,mod$scores$xscores[,1:numb_cc])
  Zy=cbind(mod$data$Zy,mod$scores$yscores[,1:numb_cc])
  
  for(i in 1:numb_cc){
    mod$p_values[i]=(sum(replicate(B,perm_and_cc(X,Y,n))>=mod$cor[i])+1)/(B+1)
    if((mod$p_values[i]>= alpha_max) | (i== numb_cc)) return(mod)
    X=residualize(mod$data$X,Zx)
    Y=residualize(mod$data$Y,Zy)
  }
  
  return(mod)
}