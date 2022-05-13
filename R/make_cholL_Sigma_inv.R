make_Z = function(term,data,K){
  stopifnot(all(unique(data[[term]]) %in% rownames(K)))
  data[[term]] = factor(data[[term]],levels = rownames(K))
  m = model.matrix(formula(sprintf('~0+%s',term)),data = data)
  colnames(m) = sub(term,'',colnames(m))
  m = m[,rownames(K)]
  m
}
khatri_rao_rows = function(X,Y) {
  stopifnot(nrow(X) == nrow(Y))
  XY = X[,rep(1:ncol(X),each = ncol(Y))] * Y[,rep(1:ncol(Y),ncol(X))]
  colnames(XY) = paste(rep(colnames(Y),ncol(X)),rep(colnames(X),each = ncol(Y)),sep=':')
  XY
}

#' make_cholSigma
#'
#' Creates a full covariance matrix for a sparse matrix of data as the sum of
#' kronecker products of covariance matrix pairs and then computes the Cholesky decomposition, and finally the inverse of this upper-triangular matrix.
#' If Y is a (sparse) matrix
#' with distribution v(Y) ~ N(0,Sigma), where Sigma = \sum K_i \otimes G_i,
#' will return the Cholesky decomposition of Sigma[j,j] where j indexes the non-NA
#' elements of Y.
#'
#' @param data data.frame of observations
#' @param phenoID name of the phenotype column in \code{data}. This vector should not have any NAs (so run \code{subset(data,!is.na(data[[phenoID]]))} first).
#' @param rowID name of the column in \code{data} that would be the row-name of Y
#' @param columnID name of the column in \code{data} that would be the column-name of Y
#' @param covariances list of Row-Column pairs of covariance matrices.list(list(Row=K,Column=G),list(Row=I,Column=R)).
#' The Row matrices should have row names, and all elements of \code{data[[rowID]]} should be present in these row names.
#' If Row is missing, it is assumed to be the identity. The Column matrices should also have row names, and all elements
#' of \code{data[[columnID]]} should be present in the row names of these matrices.
#' @param sparse If TRUE and if the Cholesky decomposition of Sigma has >25% zeros, the output will be converted to a CSparseMatrix
#'
#' @return upper-triangle matrix with the inverse of the Cholesky decomposition.
#' @export
#'
#' @examples
make_cholL_Sigma_inv = function(data,phenoID,rowID,columnID,covariances,sparse=TRUE) {
  require(Matrix)
  # recover()
  if(!phenoID %in% colnames(data)) stop(sprintf('phenoID column "%s" not found in data',phenoID))
  if(!rowID %in% colnames(data)) stop(sprintf('rowID column "%s" not found in data',rowID))
  if(!columnID %in% colnames(data)) stop(sprintf('columnID column "%s" not found in data',columnID))
  if(any(is.na(data[[phenoID]]))) stop(sprintf('NAs detected in data[["%s"]]',phenoID))
  # data = subset(data,!is.na(data[[phenoID]]))

  for(i in 1:length(covariances)) {
    cov_pair = covariances[[i]]
    if(any(names(cov_pair) %in% c('Row','Column') == F)) warning(sprintf('Possible wrong name in covariances %d',i))
    if(!is.null(cov_pair$Row)) {
      if(length(rownames(cov_pair$Row)) == 0) stop(sprintf('Row covariance %d lacks row names',i))
      if(any(data[[rowID]] %in% rownames(cov_pair$Row) == F)) stop(sprintf('rowIDs missing in Row covariance %d',i))
    }
    if(!is.null(cov_pair$Column)) {
      if(length(rownames(cov_pair$Column)) == 0) stop(sprintf('Column covariance %d lacks row names',i))
      if(any(data[[columnID]] %in% rownames(cov_pair$Column) == F)) stop(sprintf('rowIDs missing in Column covariance %d',i))
    }
  }

  Sigma = matrix(0,nrow(data),nrow(data))
  for(i in 1:length(covariances)) {
    cov_pair = covariances[[i]]
    if(is.null(cov_pair$Row)) {
      cov_pair$Row = diag(1,length(unique(data[[rowID]])))
      rownames(cov_pair$Row) = colnames(cov_pair$Row) = unique(data[[rowID]])
    }
    if(is.null(cov_pair$Column)) {
      cov_pair$Column = diag(1,length(unique(data[[columnID]])))
      rownames(cov_pair$Column) = colnames(cov_pair$Column) = unique(data[[columnID]])
    }
    Z_row = make_Z(rowID,data,cov_pair$Row)
    Z_column = make_Z(columnID,data,cov_pair$Column)
    Z = Matrix::t(Matrix::KhatriRao(t(Z_column),t(Z_row), make.dimnames=TRUE))
    KoG_Zt = khatri_rao_rows(Z_column %*% cov_pair$Column,Z_row %*% cov_pair$Row)
    stopifnot(all(colnames(Z) == colnames(KoG_Zt)))
    Sigma_i = Z %**% t(KoG_Zt)
    Sigma = Sigma + Sigma_i
    # if(mean(Sigma_i == 0) > .25 && (is(Sigma,'sparseMatrix') || all(Sigma==0))) {
    #   Sigma = Sigma + as(Sigma_i,'dgCMatrix')
    # } else{
    #   if(is(Sigma,'sparseMatrix')) Sigma = as.matrix(Sigma)
    #   Sigma = Sigma + as.matrix(Sigma_i)
    # }
  }

  if(sparse && mean(zapsmall(Sigma)==0 ) > 0.25) Sigma = as(Sigma,'dgCMatrix')
  chol_Sigma = chol(Sigma)

  chol_Sigma_inv = backsolve(chol_Sigma,diag(nrow(chol_Sigma)))

  if(sparse) {
    if(!is(chol_Sigma_inv,'sparseMatrix')) chol_Sigma_inv = as(chol_Sigma_inv,'dgCMatrix')
  #   if(is(chol_Sigma_inv,'sparseMatrix') & mean(drop0(chol_Sigma_inv,tol=1e-8)==0 ) < 0.75) {
  #     chol_Sigma_inv = as.matrix(chol_Sigma_inv)
  #     return(chol_Sigma_inv)
  #   }
  #   if(!is(chol_Sigma_inv,'sparseMatrix') & mean(drop0(chol_Sigma_inv,tol=1e-8)==0 ) > 0.75) {
  #     chol_Sigma_inv = as(chol_Sigma_inv,'dgCMatrix')
  #     return(chol_Sigma_inv)
  #   }
  # } else {
  #   if(is(chol_Sigma_inv,'sparseMatrix')) {
  #     chol_Sigma_inv = as.matrix(chol_Sigma_inv)
  #   }
  } else{
    if(is(chol_Sigma_inv,'sparseMatrix')) chol_Sigma_inv = as.matrix(chol_Sigma_inv)
  }
  return(Matrix::t(chol_Sigma_inv))
#
#   if(is(Sigma,'sparseMatrix')) {
#     chol_Sigma = chol(Sigma)
#     if(mean(chol_Sigma==0) < .75) {
#       chol_Sigma = as.matrix(chol_Sigma)
#     } else{
#       chol_Sigma = as(chol_Sigma,'dgCMatrix')
#     }
#   } else{
#     chol_Sigma = chol_c(Sigma)
#   }
#   return(chol_Sigma)
}
