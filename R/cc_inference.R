
#'cc_inference
#'@description compute p-values
#'@param mod an \code{acca} object.
#'@param B (\code{100} by default) number of random sampling
#'@param alpha (\code{0.5} by default).
#'@param numb_cc (\code{=NULL} by default).
#'@param resamp_type \code{"sign-flip"} (by default) or \code{"permutation"} 
#'@export
#'
cc_inference <- function(mod,B=100,alpha=.5,numb_cc=NULL,resamp_type="sign-flip"){
  n=nrow(mod$data$X)
  resamp_type=match.arg(resamp_type,c("sign-flip","permutation"))
  if(is.null(numb_cc)) numb_cc=length(mod$cor)
    

  
####################
  if(resamp_type=="sign-flip"){
    .permute <- function(X,Qx,nred){
      t(Qy*(1-2*rbinom(nred,1,.5))) %*%Qy%*%X
    }
  } else if(resamp_type=="permutation"){
    .permute <- function(X,Qx,nred){
      t(Qy[sample(nred),]) %*%Qy%*%X
    }
  }
   

  perm_and_cc=function(X,Y,Qx,Qy,nredx,nredy){
    ccp=.cc_core(.permute(X,Qx,nredx), .permute(Y,Qy,nredy),numb_cc = 0)
    ccp$cor[1]
  }
##################
  
  mod$p_values=rep(1,length(mod$cor))
  attr(mod$p_values,which = "B")=B
  Zx=mod$data$Zx
  Zy=mod$data$Zy
  if(is.null(Zx)) Qx=diag(nrow(mod$data$X)) else
    Qx=residualizing_matrix(Zx)$Q
  if(is.null(Zy)) Qy=diag(nrow(mod$data$Y)) else
    Qy=residualizing_matrix(Zy)$Q
  
  nredx=nrow(Qx)
  nredy=nrow(Qy)
  
  n_nuisx=min(0,ncol(Zx))
  n_nuisy=min(0,ncol(Zy))
  Zx=cbind(Zx,mod$scores$xscores[,1:numb_cc])
  Zy=cbind(Zy,mod$scores$yscores[,1:numb_cc])
  
  for(i in 1:numb_cc){
    mod$p_values[i]=(sum(replicate(B,perm_and_cc(X,Y,Qx,Qy,nredx,nredy))>=mod$cor[i])+1)/(B+1)
    if(mod$p_values[i]>=alpha) return(mod)
    Qx=residualizing_matrix(Zx[,1:n_nuisx+i])$Q
    Qy=residualizing_matrix(Zy[,1:n_nuisy+i])$Q
  }
  
  return(mod)
}