fastEMMAX_ANOVA_matrix = function(formula,data,markers,svd_matrices,mc.cores = RcppParallel::defaultNumThreads()-1,verbose = F) {
  # svd_matrices is a list with two components:
  #   svd_GR = simultaneous_diagonalize(R,Ginvsq)
  #   svd_K = svd(K)
  #   We have Sigma = kronecker(K,G) + kronecker(I,R)
  #   If we pre-multiply by kronecker(t(svd_K$u),t(svd_GR$S)), the result is a diagonal matrix
  #   We can de-correlate Y by t(svd_K$u) %*% Y %*% svd_GR$S
  #   We can de-correlate rows of X by t(svd_K$u) %*% X
  #   How to de-correlated columns of B? B %*% svd_GR$S  How to then test the columns of B?

}
