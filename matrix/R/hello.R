# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'


multi_matrix=function(A,B){
  ncol.A=ncol(A)
  ncol.B=ncol(B)
  nrow.A=nrow(A)
  #nrow.B=nrow(B)
  result=matrix(0,nrow.A,ncol.B)
  for (i in 1:nrow.A){
    for (j in 1:ncol.B){
      res.sum=0
      for(k in 1:ncol.A){
        res.sum =res.sum + A[i,k]*B[k,j]
      }
      result[i,j]=res.sum
    }
  }
  return(result)
}

multi_MMV=function(M1,M2,V,arg=1){
  V=matrix(V)
  if(arg==1){
    result=multi_matrix(M1,M2)
    result=multi_matrix(result,V)
    return(result)
  }
  else{
    result=multi_matrix(M2,V)
    result=multi_matrix(M1,result)
    return(result)
  }
}

