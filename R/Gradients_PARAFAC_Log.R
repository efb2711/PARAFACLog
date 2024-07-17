is_binary_matrix <- function(matrix) {
  return(all(matrix == 0 | matrix == 1))
}

sigmoide<-function(z){
  (1/(1+exp(-1*z)))
}

L1 <- function(A, B1, C, D1) {
  D1_pred = A %*% t(khatri_rao(C, B1))
  sum((D1 - D1_pred)^2)
}

L2 <- function(A, B2, D2) {
  M2 = sigmoide(cbind(1, A) %*% t(B2))
  sum(-D2 * log(M2) - (1 - D2) * log(1 - M2), na.rm = TRUE)
}

Loss <- function(A, B1, C, B2, D1, D2, beta) {
  L1_val = L1(A, B1, C, D1)
  L2_val = L2(A, B2, D2)
  (1-beta) * L1_val + beta * L2_val
}

LLogRegARec <- function(par, D1, D2, A, B1, C, B2, lambda, beta, s) { # Cost to estimate A
  A[,s]=par
  M2=sigmoide(cbind(1,A) %*% t(B2))
  L2=sum(-D2*log(M2)-(1-D2)*log(1-M2), na.rm = TRUE)
  D1_pred = A %*% t(khatri_rao(C,B1))
  L1 = sum((D1 - D1_pred) ^ 2)
  penalty <- lambda * sum(A^2)
  L = ((1-beta)*L1+ beta*L2)/2 + penalty/2
  return(L)
}

LLogRegB1Rec <- function(par, D1, D2, A, B1, C, B2, lambda, beta, s) { # Cost to estimate B1
  B1[,s]=par
  M2=sigmoide(cbind(1,A) %*% t(B2))
  L2=sum(-D2*log(M2)-(1-D2)*log(1-M2), na.rm = TRUE)
  D1_pred = A %*% t(khatri_rao(C,B1))
  L1 = sum((D1 - D1_pred) ^ 2)
  penalty <- lambda * sum(B1^2)
  L = ((1-beta)*L1+ beta*L2)/2 + penalty/2
  return(L)
}

LLogRegCRec <- function(par, D1, D2, A, B1, C, B2, lambda, beta, s) { # Cost to estimate C
  C[,s] = par
  M2=sigmoide(cbind(1,A) %*% t(B2))
  L2=sum(-D2*log(M2)-(1-D2)*log(1-M2), na.rm = TRUE)
  D1_pred = A %*% t(khatri_rao(C,B1))
  L1 = sum((D1 - D1_pred) ^ 2)
  penalty <- lambda * sum(C^2)
  L = ((1-beta)*L1+ beta*L2)/2 + penalty/2
  return(L)
}

LLogRegB2Rec <- function(par, D1, D2, A, B1, C, B2, lambda, beta, s) { # Cost to estimate A
  B2[,s+1]=par
  M2=sigmoide(cbind(1,A) %*% t(B2))
  L2=sum(-D2*log(M2)-(1-D2)*log(1-M2), na.rm = TRUE)
  D1_pred = A %*% t(khatri_rao(C,B1))
  L1 = sum((D1 - D1_pred) ^ 2)
  penalty <- lambda * sum(B2[,-1]^2)
  L = (1-beta)*L1+ beta*L2 + penalty
  return(L)
}

grLogRegARec <- function (par, D1, D2, A, B1, C, B2, lambda, beta, s)  { ## Gradient to estimate A
  A[,s]=par
  D1_pred = A %*% t(khatri_rao(C,B1))
  grad_L1_A = as.matrix(D1_pred - D1) %*% khatri_rao(C,B1)
  M2=sigmoide(cbind(1,A) %*% t(B2))
  grad_L2_A = (M2 - D2)%*%B2[,-1]
  grad_penalty_A <- lambda * A
  gradA = (1-beta) * grad_L1_A + beta*grad_L2_A + grad_penalty_A
  #gradA = grad_L1_A + grad_L2_A + grad_penalty_A
  grad=c(c(gradA[,s]))
  return(grad)
}

grLogRegB1Rec <- function(par, D1, D2, A, B1, C, B2, lambda, beta, s)  { ## Gradient to estimate A
  B1[,s]=par
  I = dim(A)[1]
  J = dim(B1)[1]
  K = dim(C)[1]
  D1_pred = B1 %*% t(khatri_rao(A,C))
  D1b = permnew(D1,I,J,K)
  grad_L1_B1 = as.matrix(D1_pred - D1b) %*% khatri_rao(A,C)
  grad_penalty_B1 <- lambda * B1
  gradB1 = (1-beta)*grad_L1_B1 + grad_penalty_B1
  #gradB1 = grad_L1_B1 + grad_penalty_B1
  grad=c(c(gradB1[,s]))
  return(grad)
}

grLogRegCRec <- function(par, D1, D2, A, B1, C, B2, lambda, beta, s) { ## Gradient to estimate C
  C[,s]=par
  I = dim(A)[1]
  J = dim(B1)[1]
  K = dim(C)[1]
  D1_pred = C %*% t(khatri_rao(B1,A))
  D1c = permnew(D1,J,K,I)
  grad_L1_C = as.matrix(D1_pred - D1c) %*% khatri_rao(B1,A)
  grad_penalty_C <- lambda * C
  gradC = (1-beta)*grad_L1_C + grad_penalty_C
  #gradC = grad_L1_C + grad_penalty_C
  grad=c(c(gradC[,s]))
  return(grad)
}

grLogRegB2Rec <- function (par, D1, D2, A, B1, C, B2, lambda, beta, s)  { ## Gradient to estimate B2
  B2[,s+1]=par
  M2=sigmoide(cbind(1,A) %*% t(B2))
  grad_L2_B2 = t(M2 - D2)%*%A
  grad_penalty_B2 <- lambda * B2[,-1]
  gradB2 = beta*grad_L2_B2 + grad_penalty_B2
  #gradB2 = grad_L2_B2 + grad_penalty_B2
  grad=c(c(gradB2[,s]))
  return(grad)
}
