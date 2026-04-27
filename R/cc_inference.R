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
cc_inference <-  function(mod,B=1000,stat_test="Roy",alpha_max=.5,numb_cc=NULL,resamp_type="sign-flip",light=FALSE,test_type="resampling"){
  mod$call$cc_inference=match.call()
  n=nrow(mod$data$X)
  test_type=match.arg(test_type,c("resampling","parametric"))
  if(test_type=="resampling") {
    stat_test=match.arg(stat_test, c("Roy","Wilks")) 
    resamp_type=match.arg(resamp_type, c("rotation","sign-flip","permutation"))
  } else if(test_type=="parametric") {
    stat_test="Wilks"
    light = FALSE
  }
  if(is.null(numb_cc)) numb_cc=length(mod$cor)
  
  if(test_type=="resampling") {
    if(light){
      mod=.cc_inference_residuals(mod,B,stat_test,alpha_max,numb_cc,resamp_type)
    }
    else if(!light){
      mod=.cc_inference_orthogonal(mod,B,stat_test,alpha_max,numb_cc,resamp_type)
    }
  } 
  else if(test_type=="parametric"){
      mod=.cc_inference_parametric(mod,alpha_max,numb_cc)
    }
  return(mod)
}

####################

.cc_inference_orthogonal <- function(mod,B,stat_test,alpha_max,numb_cc,resamp_type){
  resid_matrix <- function(Zy)
    residualizing_matrix(Zy)$Q
  
  if(resamp_type=="rotation") {
    .permute <- function(X, Qx, nred) {
      M = matrix(rnorm(nred*nred),nrow=nred,ncol=nred)
      O = qr.Q(qr(M,LAPACK = FALSE))
      flipsign = which(rbinom(nred, 1, 0.5) == 1)
      if(length(flipsign) > 0) O[, flipsign] = -O[, flipsign]
      Qx%*%O%*%X
    }
  } else if(resamp_type=="sign-flip") {
    .permute <- function(X,Qx,nred) {
      (1-2*rbinom(nred,1,.5))*Qx%*%X # X=t(Qx)%*%mod$data$X
    }
  } else if(resamp_type=="permutation") {
    .permute <- function(X,Qx,nred) {
      Qx[,sample(nred)]%*%X # X=t(Qx)%*%mod$data$X
    }
  }
  
  perm_and_cc <- function(X,Y,Qx,Qy,nredx,nredy,stat_test){
    ccp=.cc_core(.permute(X,Qx,nredx),Y,numb_cc = 0)
    if (stat_test=="Roy") {
      return(ccp$cor[1])
    } else if (stat_test=="Wilks") {
      return(-sum(log(1-ccp$cor^2)))
    }
  }
  
  if(is.null(mod$data$Zx)) {
    Qx=diag(nrow(mod$data$X))
    X=mod$data$X
  } else {
    Qx=resid_matrix(mod$data$Zx)
    X=t(Qx)%*%mod$data$X
  }
  
  if(is.null(mod$data$Zy)) {
    Qy=diag(nrow(mod$data$Y))
    Y=mod$data$Y
  } else {
    Qy=resid_matrix(mod$data$Zy)
    Y=t(Qy)%*%mod$data$Y
  }
  
  nredx=ncol(Qx)
  nredy=ncol(Qy)
  
  mod$p_values=rep(1,length(mod$cor))
  attr(mod$p_values,which="B")=B
  
  n_nuisx = if(is.null(mod$data$Zx)) 0 else ncol(mod$data$Zx)
  n_nuisy = if(is.null(mod$data$Zy)) 0 else ncol(mod$data$Zy)
  Zx=cbind(mod$data$Zx,mod$scores$xscores[,1:numb_cc])
  Zy=cbind(mod$data$Zy,mod$scores$yscores[,1:numb_cc])
  
  for(i in 1:numb_cc) {
    if (stat_test=="Roy") {
      obs_stat=mod$cor[i]
    } else if (stat_test=="Wilks") {
      obs_stat=-sum(log(1-mod$cor[i:length(mod$cor)]^2))
    }
    perm_stats=replicate(B,perm_and_cc(X,Y,Qx,Qy,nredx,nredy,stat_test))
    mod$p_values[i]=(sum(perm_stats >= obs_stat)+1)/(B+1)
    if(mod$p_values[i] >= alpha_max) return(mod)
    
    Qx=resid_matrix(Zx[,1:(n_nuisx+i),drop=FALSE])
    Qy=resid_matrix(Zy[,1:(n_nuisy+i),drop=FALSE])
    nredx = ncol(Qx)
    nredy = ncol(Qy)
    X=t(Qx)%*%mod$data$X
    Y=t(Qy)%*%mod$data$Y
  }
  
  return(mod)
}

####################

.cc_inference_residuals <- function(mod,B,stat_test,alpha_max,numb_cc,resamp_type){
  n=nrow(mod$data$X)
  
  if(resamp_type=="rotation") {
    .permute <- function(X,n) {
      M=matrix(rnorm(n*n),nrow=n,ncol=n)
      O = qr.Q(qr(M, LAPACK = FALSE))
      flipsign = which(rbinom(n, 1, 0.5) == 1)
      if(length(flipsign) > 0) O[, flipsign] = -O[, flipsign]
      O%*%X
    }
  } else if(resamp_type=="sign-flip") {
    .permute <- function(X,n) {
      (1-2*rbinom(n,1,.5))*X # X=IH%*%mod$data$X
    }
  } else if(resamp_type=="permutation") {
    .permute <- function(X,n) {
      X[sample(n),] # X=IH%*%mod$data$X
    }
  } 
  
  perm_and_cc <- function(X,Y,n,stat_test) {
    ccp=.cc_core(.permute(X,n),Y,numb_cc=0)
    if (stat_test=="Roy") {
      return(ccp$cor[1]) 
    } else if (stat_test=="Wilks") {
      return(-sum(log(1-ccp$cor^2)))
    }
  }
  
  if(is.null(mod$data$Zx)) {
    X=mod$data$X
  } else {
    X=residualize(mod$data$X,mod$data$Zx)
  }
  
  if(is.null(mod$data$Zy)) {
    Y=mod$data$Y
  } else {
    Y=residualize(mod$data$Y,mod$data$Zy)
  }
  
  mod$p_values=rep(1,length(mod$cor))
  attr(mod$p_values,which="B")=B
  
  n_nuisx = if(is.null(mod$data$Zx)) 0 else ncol(mod$data$Zx)
  n_nuisy = if(is.null(mod$data$Zy)) 0 else ncol(mod$data$Zy)
  Zx=cbind(mod$data$Zx,mod$scores$xscores[,1:numb_cc])
  Zy=cbind(mod$data$Zy,mod$scores$yscores[,1:numb_cc])
  
  for(i in 1:numb_cc) {
    if (stat_test=="Roy") {
      obs_stat=mod$cor[i]
    } else if (stat_test=="Wilks") {
      obs_stat=-sum(log(1-mod$cor[i:length(mod$cor)]^2))
    }
    perm_stats <- replicate(B, perm_and_cc(X,Y,n,stat_test))
    mod$p_values[i]=(sum(perm_stats >= obs_stat)+1)/(B+1)
    if(mod$p_values[i] >= alpha_max) return(mod)
    
    X = residualize(mod$data$X,Zx[,1:(n_nuisx+i),drop=FALSE])
    Y = residualize(mod$data$Y,Zy[,1:(n_nuisy+i),drop=FALSE])
    
  }
  return(mod)
}

####################

.cc_inference_parametric <- function(mod,alpha_max,numb_cc){
  
  n=nrow(mod$data$X)
  r = if(is.null(mod$data$Zx)) 0 else qr(mod$data$Zx)$rank
  n_eff=n-r
  p=ncol(mod$data$Y)
  q=ncol(mod$data$X)
  bartlett=-(n_eff-1-1/2*(p+q+1))
  
  mod$p_values=rep(1,length(mod$cor))
  attr(mod$p_values,which ="test_type")="parametric"
  
  for (i in 1:numb_cc) {
    lambda=log(prod(1-(mod$cor[i:length(mod$cor)])^2))
    df=(p-i+1)*(q-i+1)
    mod$p_values[i]=pchisq(bartlett*lambda,df,lower.tail=FALSE)
    if(mod$p_values[i] >= alpha_max) return(mod)
  }
  return(mod)
}

####################
