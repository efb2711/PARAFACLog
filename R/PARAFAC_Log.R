
# PARAFAC_Logit
#'
#' This function fits a PARAFAC-Log model to a coupled data set consisting of
#' a three-way array of real data and a binary data matrix that share a mode.
#'
#' @param X The first mode unfolding of a three-way array of real-valued data.
#' @param Y A binary data matrix.
#' @param I Number of elements in the first mode in X.
#' @param J Number of elements in the second mode in X.
#' @param K Number of elements in the third mode in X.
#' @param L Number of columns in Y.
#' @param inames: Optional, names for the elements in the first mode of X.
#' @param xnames: Optional, names for the elements in the second/third mode of X.
#' @param ynames: Optional, names of the variables of Y.
#' @param R Number of components in the PARAFAC model.
#' @param alfa Weighting parameter.
#' @param tolerance Convergence tolerance.
#' @param maxiter Maximum number of iterations.
#' @param penalization Learning rate,
#' @param OptimMethod Optimization method to use (default is "CG").
#
#' @return A list with the following elements:
#' \item{A}{Component matrix A.}
#' \item{B1}{Component matrix B1.}
#' \item{C}{Component matrix C.}
#' \item{B2}{Factor matrix B2.}
#' \item{b0}{Intercept.}
#' \item{L}{Final value of the loss function.}
#' \item{Expected}{Matrix of expected values.}
#' \item{Predictions}{Matrix of binary predictions.}
#'
#'@examples
#' Example usage:
#' X <- array(rnorm(1000), dim = c(10, 10, 10))
#' Y <- matrix(rbinom(100, 1, 0.5), nrow = 10, ncol = 10)
#' result <- PARAFAC_Logit_Fit(X, Y, 10, 10, 10, 10)
#' @export
PARAFAC_Logit <- function(X, Y, I, J, K, L, inames=NULL, xnames=NULL, ynames=NULL,
                          R=2, alfa = 0.5, tolerance=0.000005, maxiter=100,
                          penalization=0.1, OptimMethod="CG")
{
  # INPUT
  #    X (I x J x K): 3D-data matrix (represented as a 3D-array) [ X is concatenated version of X3D ]
  #    Y (I x L): 2D-data matrix
  #    I: number of elements of first mode of 3D/2D (the common mode: rows)
  #    J: number of elements of second mode of 3D (columns 3D)
  #    K: number of elements of third mode of 3D (slabs)
  #    L: number of elements of second mode of 2D (columns 2D)
  #    r: rank of the model PARAFAC
  #    tolerance: tolerance value
  #    alfa: alfa value
  #    maxiter: number of runs
  X = as.matrix(X)
  if( R > min(I,J,K,L) ) stop ("rank should be an integer between 1 and " , min(I,J,K,L))

  if ( (alfa < 0) || (alfa >= 1) ) stop ("Alfa should be between 0 and 1 (but not 0 or 1)")

  if (!is_binary_matrix(Y)) stop("La matriz de entrada no es binaria.")

  dimnames=paste("Comp.", 1:R)
  if (is.null(inames)){
    inames = paste("i",1:I,sep="")
  }
  if (is.null(xnames)){
    for (k in 1:K)
      for (j in 1:J){
        xnames=c(xnames, paste(k, j,sep="."))
      }
  }
  if (is.null(ynames)){
    ynames = paste("Y",1:L,sep="")
  }
  rownames(X) = inames
  colnames(X) = xnames
  #ynames= c("Fear to be refused","Kindness","Importance of others’ judgments","Altruism","Neuroticism","Being strict to oneself","Low selfesteem","Conscientiousness","Depression")
  result=list()
  result$Method="PARAFAC-Logit"
  result$X=X
  result$Y=Y
  result$tolerance=tolerance
  result$maxiter=maxiter
  result$penalization=penalization
  myfit=PARAFAC_Logit_Fit(X=X, Y=Y,I=I, J=J, K=K, L=L, R=R, alfa = alfa, tolerance=tolerance, maxiter=maxiter,
                    penalization=penalization,OptimMethod=OptimMethod)
  result$X = X
  result$Y = Y
  result$A = myfit$A
  result$B1 = myfit$B1
  result$C = myfit$C
  result$B2 = myfit$B2
  result$b0 = myfit$b0
  result$Expected = myfit$Expected
  result$Predictions=myfit$Predictions
  colnames(result$Predictions)=ynames
  class(result)="PARLog"
  return(result)
}
