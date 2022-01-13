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

#' multiply two matrix
#' @param M1 a matrix
#' @param M2 a matrix
#' @return Multiply M1, M2
#' @examples M1=matrix(c(1,2,3,9),ncol=2)
#' M2=matrix(c(9,2,3,9),ncol=2)
#' multi_matrix(M1,M2)
multi_matrix=function(M1,M2){
  ncol.M1=ncol(M1)
  ncol.M2=ncol(M2)
  nrow.M1=nrow(M1)
  #nrow.M2=nrow(M2)
  result=matrix(0,nrow.M1,ncol.M2)
  for (i in 1:nrow.M1){
    for (j in 1:ncol.M2){
      res.sum=0
      for(k in 1:ncol.M1){
        res.sum =res.sum + M1[i,k]*M2[k,j]
      }
      result[i,j]=res.sum
    }
  }
  return(result)
}


#' multiply two matrix and a vector
#' @param M1 a matrix
#' @param M2 a matrix
#' @param V a vector
#' @param arg The order of the multiplication, 1 first multiply  M1 M2, 2 first multiply M2 V
#' @return Multiply M1, M2, V
#' @examples M1=matrix(c(1,2,3,9),ncol=2)
#' M2=matrix(c(9,2,3,9),ncol=2)
#' V=c(1,2)
#' multi_MMV(M1,M2,V,arg=1)
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

