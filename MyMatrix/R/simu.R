


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
