#' estimate parameters for the logisitic regression
#' @param Y a vector with n 0-1 entry
#' @param X data
#' @param mod estimate mod
#' @param alpha alpha for the estimate
#' @return estimated parameters, CI, log Likelihood
op.log=function(Y,X,mod,alpha=1){
  X.nrow=nrow(X)
  X=cbind(matrix(rep(1,X.nrow)),as.matrix(X))
  X.ncol=ncol(X)
  v.beta=matrix(rep(0,X.ncol))
  error=10
  l=c(0)
  p=1/(1+exp(-X%*%v.beta))
  while(error>0.1){
    beta.old=v.beta
    if (mod=='Newton'){
      #W=-p%*%t(1-p)
      W=diag(c(p*(1-p)))
      At=-t(X)%*%W%*%X
      alpha=1
    }
    else if(mod=='Gradient'){
      At=diag(rep(1,X.ncol))
      colnames(At)=colnames(X)
    }
    v.beta=v.beta-alpha*solve(At)%*%t(X)%*%(Y-p)
    p=1/(1+exp(-X%*%v.beta))
    l=append(l,sum(Y*log(min(p,1-1e-100)))+sum((1-Y)*log(max(1-p,1e-100))))
    #l=append(l,sum(Y*log(p))+sum((1-Y)*log(1-p)))
    error=abs(l[length(l)]-l[length(l)-1])
  }
  M=t(X)%*%diag(c(p*(1-p)))%*%X
  result=data.frame(value=v.beta,
                    CI.L=v.beta+qnorm(0.025,0,1)*matrix(sqrt(1/diag(M))),
                    CI.H=v.beta+qnorm(0.975,0,1)*matrix(sqrt(1/diag(M))))
  return(list(result,l[-1]))
}


#' estimate parameters of EM algorithm
#' @param data a vector of data
#' @param init.value initial value for the vector p.
#' @return estimated parameters of p
#' @example EM(c(6,4,55,35),c(.1,.1,.8))
EM=function(data,init.value){
  na=data[1]
  nab=data[2]
  nb=data[3]
  no=data[4]
  n=na+nab+nb+no
  pa=init.value[1]
  pb=init.value[2]
  po=init.value[3]
  eps=1
  p.old=0
  p=1
  while( eps>0.001 & abs(p-p.old)>0.001){
    pa.old=pa
    pb.old=pb
    po.old=po
    p.old=p
    maa=na*(pa^2)/(pa^2+2*pa*po)
    mao=na*(2*pa*po)/(pa^2+2*pa*po)
    mbb=nb*(pb*pb)/(pb^2+2*pb*po)
    mbo=nb*(2*pb*po)/(pb^2+2*pb*po)
    mab=nab
    moo=no
    pa=(2*maa+mao+mab)/2/n
    pb=(2*mbb+mbo+mab)/2/n
    po=(2*moo+mao+mbo)/2/n
    eps=sqrt((pa-pa.old)^2+(pb-pb.old)^2+(po-po.old)^2)
    p=(pa^2+2*pa*po)^(na)*(pb^2+2*pb*po)^(nb)*(po^2)^(no)
  }
  return(list(c(pa=pa,pb=pb,po=po,p=p)))
}
