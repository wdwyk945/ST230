


#' generate MVN samples
#' @param mu a n dimensional numeric vector
#' @param sigma a n*n positive definite matrix
#' @param N a numeric value
#' @return Sample data from MVN(mu,sigma) for N times
#' @examples sigma=matrix(c(1,.5,.5,.5,.5,1,5,.5,.5,.5,1,.5,.5,.5,.5,1),ncol=4)
#' data=simu(c(1,2,3,4), sigma,100)
#' summary(data)
#' cov(data)
simu=function(mu, sigma,N){
  data=data.frame(matrix(rep(0,N*length(mu)),nrow=N))

  for(i in 1:N){
    Z=rnorm(length(mu),0,1)
    X=mu+Z%*%chol(sigma)
    data[i,]=X
  }
  return(data)
}

#' Cholesky decomposition
#' @param sigma a n*n positive definite matrix
#' @return Cholesky decomposition of sigma
#' @examples sigma=matrix(c(1,.5,.5,.5,.5,1,5,.5,.5,.5,1,.5,.5,.5,.5,1),ncol=4)
#' chol(sigma)
chol=function(sigma){
  ncol.sigma=ncol(sigma)
  nrow.sigma=nrow(sigma)
  L=matrix(rep(0,
               ncol.sigma*nrow.sigma),
           ncol=ncol.sigma)
  l11=sqrt(sigma[1,1])
  l=sigma[1,2:ncol.sigma]/l11
  M=sigma[2:nrow.sigma,2:ncol.sigma]-l%*%t(l)
  L[1,1]=l11
  L[1,2:ncol.sigma]=l
  L[2:nrow.sigma,1]=0
  if(nrow.sigma==2){
    L[2,2]=sqrt(M)
  }
  else{
    L[2:nrow.sigma,2:ncol.sigma]=chol(M)
  }
  return(L)
}

#QR.decom=function(X){
#  X.ncol=ncol(X)
#  X.nrow=nrow(X)
#  U=matrix(rep(0,X.ncol*X.nrow),ncol=X.ncol)
#  V=matrix(rep(0,X.ncol*X.nrow),ncol=X.ncol)
#  U[1,]=X[1,]/sqrt(sum(X[1,]^2))
#  V[1,]=X[1,]
#  for (i in 2:X.nrow){
#    u=matrix(rep(0,X.ncol),ncol=X.ncol)
#  for (j in 1:(i-1)){
#    u=u+U[j,]%*%matrix(X[i,])%*%U[j,]
#  }
#  V[i,]=X[i,]-u
#  U[i,]=V[i,]/sqrt(sum(V[i,]^2))
#}
#R=matrix(rep(0,X.nrow*X.nrow),ncol=X.nrow)
#for (i in 1:X.nrow){
#  for (j in 1:X.ncol){
#    R[i,j]=U[i,]%*%X[j,]
#  }
#}
#R=U%*%t(X)
#return(list(U,V,R))
#t(result[[1]])%*%result[[3]]
#}


#' OLS by QR decomposition
#' @param Y a vector with n entry
#' @param X a n * p matrix
#' @return estimated parameters for OLS regression by QR decomposition
#' @examples Y=rnorm(100,10,1)
#' X=matrix(rnorm(10000,10,1),ncol=100)
#' QR_OLS(Y,X)
QR_OLS=function(Y,X){
  #result=QR.decom(X)
  result=qr(X)
  Q=qr.Q(result)
  R=qr.R(result)
  beta=solve(R)%*%t(Q)%*%Y
  return(beta)
}

#' OLS by SVD decomposition
#' @param Y a vector with n entry
#' @param X a n * p matrix
#' @return estimated parameters for OLS regression by SVD decomposition
#' @examples Y=rnorm(100,10,1)
#' X=matrix(rnorm(10000,10,1),ncol=100)
#' SVD_OLS(Y,X)
SVD_OLS=function(Y,X){
  result=svd(X)

  u=result$u
  sigma.inv=solve(diag(result$d))
  v=s=result$v
  beta=v%*%sigma.inv%*%t(u)%*%Y
  return(beta)
}


