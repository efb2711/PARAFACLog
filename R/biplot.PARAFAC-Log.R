#' Biplot for PARAFAC-Logit results
#'
#' This function generates a biplot from the results of a PARAFAC-Logit model,
#' including various statistics and contributions.
#'
#' @param PARLog A list containing the results of a PARAFAC-Logit model,
#' including factor matrices, predictions, and other relevant data.
#' @return A list containing the biplot information, including row and column
#' coordinates, scale factor, eigenvalues, inertia, contributions, and various
#' fit statistics.
#' @examples
#' # Example usage:
#' # result <- PARAFAC_Logit(X, Y, I, J, K, L, R=2, alfa=0.5, tolerance=0.000005, maxiter=100, penalization=0.1, OptimMethod="CG")
#' # biplot <- Biplot.PARAFACLog(result)
#' @export
#'
Biplot.PARAFACLog <- function(PARLog)
{
  X=PARLog$X
  Y=PARLog$Y

  I=dim(X)[1]
  JK=dim(X)[2]
  L=dim(Y)[2]
  R=dim(PARLog$A)[2]

  Biplot = list()
  Biplot$Title = " PARAFAC-Logit - Biplot"
  Biplot$Type = "PARAFAC-Logit"
  Biplot$Dimension=R
  Biplot$alpha=0
  #Biplot$Initial_Transformation=plsr$Initial_Transformation
  Biplot$ncols=JK
  Biplot$nrows=I
  Biplot$dim=R

  a=PARLog$A
  b=khatri_rao(PARLog$C, PARLog$B1)
  sca = sum(a^2)
  scb = sum(b^2)
  sca = sca/I
  scb = scb/JK
  scf = sqrt(sqrt(scb/sca))
  a = a * scf
  b = b/scf

  Biplot$RowCoordinates = a
  Biplot$ColCoordinates = b
  Biplot$Scale_Factor = scf
  BC = khatri_rao(PARLog$C, PARLog$B1)
  Cont=CalculateContributions(PARLog$X, PARLog$A,BC )
  Biplot$EigenValues = Cont$Fit
  Biplot$Inertia = Cont$Fit * 100
  Biplot$CumInertia = cumsum(Biplot$Inertia)
  Biplot$RowContributions = Cont$RowContributions
  Biplot$ColContributions = Cont$ColContributions
  Biplot$Structure = Cont$Structure
  class(Biplot)="ContinuousBiplot"

  nulllinterm = matrix(1,I,1) %*% matrix(PARLog$b0, nrow=1)
  nullfitted = exp(nulllinterm)/(1 + exp(nulllinterm))
  nullDeviance = -2 * apply((Y * log(nullfitted) + (1 - Y) * log(1 - nullfitted)), 2,sum)
  modelDeviance = -2 * apply((Y * log( PARLog$Expected) + (1 - Y) * log(1 -  PARLog$Expected)), 2, sum)
  modelDif=(nullDeviance - modelDeviance)
  modeldf=R
  modelp=1-pchisq(modelDif, df =  modeldf)

  CoxSnell=1-exp(-1*modelDif/I)
  Nagelkerke=CoxSnell/(1-exp((nullDeviance/(-2)))^(2/I))
  MacFaden=1-(modelDeviance/nullDeviance)
  residuals=Y-PARLog$Expected


  pred2 = matrix(as.numeric(PARLog$Expected > 0.5), I, L)
  acier = matrix(as.numeric(round(Y) == pred2), I, L)
  presences=apply(Y, 2, sum)
  absences=I-presences
  sens = apply((acier==1) & (Y==1), 2, sum)/presences
  spec = apply((acier==1) & (Y==0), 2, sum)/absences
  totsens = sum((acier==1) & (Y==1))/sum(presences)
  totspec = sum((acier==1) & (Y==0))/sum(absences)


  SSRes=apply((residuals^2), 2,sum)
  SSTot=apply((Y^2), 2, sum)
  R2 = 1 - SSRes/SSTot

  Res=list()
  Res$Biplot="Binary Logistic (from PARAFAC-Log)"
  Res$Type= "Binary Logistic (from PARAFAC-Log)"
  Res$ColumnParameters=cbind(PARLog$b0, PARLog$B2)
  rownames(Res$ColumnParameters)=colnames(PARLog$Predictions)
  # #Res$ColumnParameters=Res$ColumnParameters
  Res$NullDeviances=nullDeviance
  Res$ModelDeviances=modelDeviance
  Res$ModelDevianceTotal=sum(Res$ModelDeviances)
  Res$NullDevianceTotal=sum(Res$NullDeviances)
  Res$Deviances=modelDif
  Res$Dfs=modeldf
  Res$pvalues=modelp
  Res$Bonferroni=(modelp * L)* ((modelp * L)<=1) + (((modelp * L)>1))
  Res$CoxSnell=CoxSnell
  Res$Nagelkerke=Nagelkerke
  Res$MacFaden=1-(Res$ModelDeviances/Res$NullDeviances)
  Res$R2=R2
  Prediction=matrix(as.numeric(PARLog$Expected>0.5), I,L)
  Correct=matrix(as.numeric(Y==Prediction), I,L)
  Res$Scale_Factor = scf
  Res$DevianceTotal=sum(Res$Deviances)
  nn=I*L
  Res$TotCoxSnell=1-exp(-1*Res$DevianceTotal/nn)
  Res$TotNagelkerke=Res$TotCoxSnell/(1-exp((Res$NullDevianceTotal/(-2)))^(2/nn))
  Res$TotMacFaden=1-(Res$ModelDevianceTotal/Res$NullDevianceTotal)
  Res$TotalDf=L*R
  SSResT=sum(residuals^2)
  SSTotT=sum(Y^2)
  Res$TotR2 = 1 - SSResT/SSTotT

  Res$PercentsCorrec=apply(Correct, 2, sum)/I
  Res$TotalPercent=sum(Correct)/(I*L)
  Res$Sensitivity=sens
  Res$Specificity=spec
  Res$TotalSensitivity=totsens
  Res$TotalSpecificity=totspec
  Res$p=1-pchisq(Res$DevianceTotal, df = Res$TotalDf)

  dd = sqrt(rowSums(cbind(1,Res$ColumnParameters[, 2:(R + 1)])^2))
  Res$Loadings = diag(1/dd) %*% Res$ColumnParameters[, 2:(R + 1)]
  Res$Tresholds = Res$ColumnParameters[, 1]/modelDif
  Res$Communalities = rowSums(Res$Loadings^2)
  Res$ColumnParameters[, 2:(R + 1)]=Res$ColumnParameters[, 2:(R + 1)]/scf
  Biplot$BinSupVarsBiplot=Res
  class(Biplot$BinSupVarsBiplot)="BinSupVarsBiplot"
  return(Biplot)
}

