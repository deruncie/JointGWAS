
#' Drops markers by missing rate and/or minor allele frequency
#'
#' @param mat n x m marker matrix (numeric) for n individuals and markers
#' @param max_missing markers with mor than this fraction of missing values will be dropped
#' @param min_maf markers with a minor allele frequency lower than this will be dropped
#'
#' @return matrix with columns failing either \code{max_missing} or \code{min_maf} dropped
#' @export
#'
#' @examples
drop_markers = function(mat,max_missing = 0.1,min_maf = 0.05) {
  mat[,apply(mat,2,function(x) mean(!is.na(x)) >= max_missing & (.5-abs(.5-mean(x,na.rm=T)))>min_maf)]
}
#' Fills in missing genotype values with the mean value for the marker
#'
#' Note: this is a quick work-around, but ideally with dense markers more effective imputation should be used (e.g. Beagle)
#'
#' @param mat n x m marker matrix (numeric) for n individuals and markers
#'
#' @return input matrix with NAs replaced with colMeans.
#' @export
#'
#' @examples
fill_missing = function(mat) {
  means = colMeans(mat,na.rm=T)
  mat[is.na(mat)] = matrix(means,nrow = nrow(mat),ncol = ncol(mat),byrow=T)[is.na(mat)]
  mat
}


`%**%` = function(X1,X2){
  if(is.null(X1)) return(X2)
  return(matrix_multiply_toDense(X1,X2))
  # if(inherits(X1,'dgCMatrix') && inherits(X2,'matrix')) return(SxD(X1,X2))
  # if(inherits(X1,'dgCMatrix') && inherits(X2,'dgCMatrix')) return(SxS(X1,X2))
  # if(inherits(X1,'matrix') && inherits(X2,'matrix')) return(X1 %*% X2)
  # if(inherits(X1,'data.frame') || inherits(X2,'data.frame')) return(as.matrix(X1) %**% as.matrix(X2))
  # return(as.matrix(X1 %*% X2))
}


simultaneous_diagonalize = function(A,B) {
  sB = svd(B)
  Binvsq = t(1/sqrt(sB$d)*t(sB$u))
  sBAB = svd(t(Binvsq) %*% A %*% Binvsq)
  O = sBAB$u
  S = Binvsq %*% O
  return(list(u=S,d=sBAB$d))
}


# base_design_matrix = function(formula,data) {
#   m = lm(formula,data,model=T)
#   X_base = model.matrix(m)
#   # if X_base is not full rank, we need to drop columns
#   if(m$qr$rank < ncol(X_base)) {
#     assign = attributes(X_base)$assign
#     contrasts = attributes(X_base)$contrasts
#     X_base = X_base[,m$qr$pivot[1:m$qr$rank]]
#     attributes(X_base)$assign = assign[m$qr$pivot[1:m$qr$rank]]
#     attributes(X_base)$contrasts = contrasts
#   }
#   X_base
# }

# construct_tests = function(X_base,mat,data,genotypeID) {
#   if(!genotypeID %in% colnames(data)) stop(sprintf('genotypeID %s missing from data',genotypeID))
#   if(any(rownames(mat) %in% data[[genotypeID]] == F)) stop('some genotypes missing from mat (check rownames)')
#   X_test = do.call(cbind,lapply(seq_len(ncol(mat)),function(i) X_base*mat[data[[genotypeID]],i]))
# }
