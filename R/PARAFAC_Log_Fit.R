#' PARAFAC-Log Model Fit
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
#' @param R Number of components in the PARAFAC model.
#' @param alfa Weighting parameter.
#' @param tolerance Convergence tolerance.
#' @param maxiter Maximum number of iterations.
#' @param penalization Learning rate,
#' @param OptimMethod Optimization method to use (default is "CG").
#'
#' @return A list with the following elements:
#' \item{A}{Component matrix A.}
#' \item{B1}{Component matrix B1.}
#' \item{C}{Component matrix C.}
#' \item{B2}{Factor matrix B2.}
#' \item{b0}{Intercept.}
#' \item{L}{Final value of the loss function.}
#' \item{Expected}{Matrix of expected values.}
#' \item{Predictions}{Matrix of binary predictions.}
#' @examples
#' Example usage:
#' result <- PARAFAC_Logit_Fit(X=X, Y=Y, I=8, J=7, K=6, L=9, R=2)
#' @export
PARAFAC_Logit_Fit <- function(X, Y, I, J, K, L, R=2, alfa = 0.5, tolerance=0.000001, maxiter=100,
                        penalization=0.1, OptimMethod="CG")
{
    set.seed(0)
    #Calcula beta
    sum_squared_D1 = sum(X^2)
    p = colMeans(Y^2)
    l0 <- sum(-Y^2 * log(p) - (1 - Y^2) * log(1 - p))
    beta = (alfa * sum_squared_D1) / (alfa * sum_squared_D1 + (1 - alfa) * l0)

    #Ajusta b0
    b0 = matrix(0, nrow = L, ncol = 1)
    for (l in 1:L) b0[l] = RidgeBinaryLogistic(y = Y[, l], matrix(1,I, 1), penalization = penalization*beta)$beta

    A=matrix(1,nrow=I, 1)
    B1=matrix(1,nrow=J, 1)
    C=matrix(1,nrow=K, 1)
    B2 = matrix(1,nrow=L,1)
    B2=b0
    X = as.matrix(X)

    for (r in 1:R){
      #
      # inicializo según una normal N(0,1) las matrices A,B1,C,B2
      # Necesito añadir una columna de 1 en A para poder multiplicarla por B2
      #
      #parA = rnorm( I , 0 , 1 )
      EIG = eigen(X %*% t(X))
      parA = EIG$vectors[, r]
      rm(EIG)

      Z = permnew(X , I , J , K)
      EIG = eigen(Z %*% t(Z))
      parB1 = EIG$vectors[, r]
      rm(EIG)

      Z = permnew(X , J , K, I)
      EIG = eigen(Z %*% t(Z))
      parC = EIG$vectors[, r]
      # parB1 = rnorm( J , 0 , 1 )
      #parC = rnorm( K , 0 , 1 )
      parB2 = rnorm( L , 0 , 1 )
      if (r == 1){
        A[,r]=parA
        B1[,r]=parB1
        C[,r]=parC
      }else
      {
        A = cbind(A,parA)
        B1 = cbind(B1,parB1)
        C = cbind(C,parC)
      }
      B2 = cbind(B2,parB2)
      Lnew = Loss(A,B1,C,B2,X,Y,beta)
      err=1
      iter=0
      while( (err > tolerance) & (iter<maxiter)){

        iter=iter+1
        Lold=Lnew
        #Update A
        resbipA <- optim(parA,fn=LLogRegARec,gr=grLogRegARec, method=OptimMethod, D1=X, D2=Y, A=A, B1=B1, B2=B2, C=C,lambda=penalization,beta = beta, s=r)
        parA=resbipA$par
        A[,r]=parA
        #Update B1
        resbipB1 <- optim(parB1,fn=LLogRegB1Rec,gr=grLogRegB1Rec, method=OptimMethod, D1=X, D2=Y, A=A, B1=B1, B2=B2, C=C,lambda=penalization,beta = beta, s=r)
        parB1=resbipB1$par
        B1[,r]=parB1
        #Update C
        resbipC <- optim(parC,fn=LLogRegCRec,gr=grLogRegCRec, method=OptimMethod, D1=X, D2=Y, A=A, B1=B1, B2=B2, C=C,lambda=penalization,beta = beta, s=r)
        parC=resbipC$par
        C[,r]=parC
        #Update B2
        resbipB2 <- optim(parB2,fn=LLogRegB2Rec,gr=grLogRegB2Rec, method=OptimMethod, D1=X, D2=Y, A=A, B1=B1, B2=B2, C=C,lambda=penalization,beta = beta, s=r)
        parB2=resbipB2$par
        B2[,r+1]=parB2

        Lnew = Loss(A,B1,C,B2,X,Y,beta)
        err=abs(Lold-Lnew)/Lnew
        cat("\n",round(iter), round(Lnew, 3), round(err,6))
      }
    }
  # Predicciones
  # D1_pred = A %*% t(khatri_rao(C, B1))
  # M2 = sigmoide(cbind(1, A) %*% t(B2))
  # BOF3D = sum((X - D1_pred)^2)                                                                             2) / ssq2D))) * 100
  # FitPercentage3D = 100 * ((ssq3D - BOF3D) / ssq3D)
  Lin= cbind(1, A) %*% t(B2)
  Expected=exp(Lin)/(1+exp(Lin))
  Pred=matrix(as.numeric(Expected>0.5), nrow=I)
  result=list(A = A, B1 = B1, C = C, B2= B2[,-1], b0 = b0, L=Lnew, Expected=Expected,
              Predictions=Pred)
}

