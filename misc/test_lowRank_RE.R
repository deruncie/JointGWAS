library(lme4)
library(lme4qtl)

n=100
d = data.frame(ID = 1:n,G = factor(rep(1:10,each = n/10)))
d$y = rnorm(length(unique(d$G)))[d$G] + rnorm(nrow(d))

m1 = lmer(y~(1|G),d)
Z = model.matrix(~0+G,d)
K = tcrossprod(Z)
rownames(K) = colnames(K) = d$ID
m2 = relmatLmer(y~(1|ID),d,relmat = list(ID = K))

Z2 = matrix(0,nrow(d),nrow(d),dimnames = list(d$ID,d$ID))
Z2[,1:ncol(Z)] = Z

# relfac <- function(mat, method.relfac = "auto", tol.relfac.svd = 1e-10, tol.relfac.evd = 1e-10)
# {
#   # ?Matrix::chol
#   # Returned value: a matrix of class Cholesky, i.e., upper triangular: R such that R'R = x.
#   # Note that another notation is equivalent x = L L', where L is a lower triangular
#   # @ http://en.wikipedia.org/wiki/Cholesky_decomposition
#   #
#   # If the substitution is Z* = Z L, then Z*' = L' Z' = R Z'
#   # @ http://www.journalofanimalscience.org/content/88/2/497.long%5BVazquez%20et%20al.,%202010%5D
#
#   # inc
#   stopifnot(requireNamespace("Matrix", quietly = TRUE))
#
#   rmat <- switch(method.relfac,
#                  "chol" = relfac.chol(mat),
#                  "svd" = relfac.svd(mat, tol.relfac.svd),
#                  "evd" = relfac.evd(mat, tol.relfac.evd),
#                  "auto" = {
#                    # filters
#                    mat.duplicated <- any(duplicated(llply(1:ncol(mat), function(i) mat[, i])))
#
#                    # auto
#                    if(mat.duplicated) {
#                      rmat <- relfac.evd(mat, tol.relfac.evd)
#                    } else {
#                      rmat <- suppressWarnings(try({
#                        relfac.chol(mat)
#                      }, silent = TRUE))
#
#                      if(class(rmat) == "try-error") {
#                        rmat <- relfac.evd(mat, tol.relfac.evd)
#                      }
#                    }
#
#                    rmat
#                  },
#                  stop())
#
#   rownames(rmat) <- rownames(mat)
#   colnames(rmat) <- colnames(mat)
#
#   return(rmat)
# }
relfac.chol <- function(mat)
{
  print('asdf')
  # recover()
  as(Matrix::t(mat),'dgCMatrix')
  # rmat = mat[,colSums(mat^2)>0,drop=FALSE]
  # rmat[,colSums(rmat^2)>0,drop=FALSE]
}
assignInNamespace('relfac.chol',relfac.chol,ns = 'lme4qtl')
# assignInNamespace('relfac',relfac,ns = 'lme4qtl')
m3 = relmatLmer(y~(1|ID),d,relmat = list(ID = Z2),method.relfac='chol')


library(lme4qtl)
d = data.frame(ID = 1:nrow(PY),y=PY[,1])
Z = matrix(0,nrow(d),nrow(d),dimnames = list(d$ID,d$ID))
Z[,1:ncol(PX_design)] = PX_design
m = relmatLmer(y~0+(1|ID),d,relmat = list(ID = Z),method.relfac='chol')


X_main = c()
X_int = c()
for(i in which(index)) {
  print(i)
  data$X = markers[data[[genotypeID]],i]
  nas = is.na(data$X)
  if(any(nas)) {
    data$X[nas] = mean(data$X,na.rm=T)
  }
  data$X = data$X - as.numeric(names(sort(table(data$X),decreasing = T)[1]))
  X_design = Matrix::sparse.model.matrix(marker_formula,data)
  assign = attr(X_design,'assign')
  X_main = cbind(X_main,X_design[,assign==1])
  X_int = cbind(X_int,X_design[,assign == 2])
}
j = Matrix::rowSums(X_design != 0) != 0

X_design = X_main
cVi_design = partial_matrix_multiply_toDense(cholL_Sigma_inv,X_design,j)
cVi_design_sparse = as(cVi_design,'dgCMatrix')
PX_design = cVi_design - cVi_Xcov_Xcovt_Vi_Xcov_inv %*% (t_cVi_Xcov %**% cVi_design_sparse)
colnames(PX_design) = colnames(X_design)
PX_main = PX_design

X_design = X_int
cVi_design = partial_matrix_multiply_toDense(cholL_Sigma_inv,X_design,j)
cVi_design_sparse = as(cVi_design,'dgCMatrix')
PX_design = cVi_design - cVi_Xcov_Xcovt_Vi_Xcov_inv %*% (t_cVi_Xcov %**% cVi_design_sparse)
colnames(PX_design) = colnames(X_design)
PX_int = PX_design

d = data.frame(ID = 1:nrow(PY),y=PY[,1])
Z = matrix(0,nrow(d),nrow(d),dimnames = list(d$ID,d$ID))
Z[,1:ncol(PX_int)] = PX_int
m = relmatLmer(y~0+PX_main+(1|ID),d,relmat = list(ID = Z),method.relfac='chol')
m = relmatLmer(y~0+(1|ID),d,relmat = list(ID = Z),method.relfac='chol')

X_design = cbind(X_main,X_int)
cVi_design = partial_matrix_multiply_toDense(cholL_Sigma_inv,X_design,j)
cVi_design_sparse = as(cVi_design,'dgCMatrix')
PX_design = cVi_design - cVi_Xcov_Xcovt_Vi_Xcov_inv %*% (t_cVi_Xcov %**% cVi_design_sparse)
colnames(PX_design) = colnames(X_design)
assign = c(rep(1,ncol(X_main)),rep(2,ncol(X_int)))
object = lm.fit(PX_design,PY)
object$df.residual = object$df.residual - b_cov
object$terms = terms(marker_formula)
object$assign = assign
class(object) = 'lm'
anova(object)
