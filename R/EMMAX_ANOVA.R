collect_results = function(results_list) {
  traits = rownames(results_list[[1]]$anova)
  results = lapply(traits,function(trait) {
    results_trait = list(
      beta_hats = do.call(dplyr::bind_rows,lapply(results_list,function(x) x$beta_hats[,trait])),
      SEs = do.call(dplyr::bind_rows,lapply(results_list,function(x) x$SEs[,trait])),
      anova = do.call(dplyr::bind_rows,lapply(results_list,function(x) data.frame(x$anova[trait,,drop=FALSE])))
    )
  })
  names(results) = traits
  results
}

anova_table = function(anova) {
  if(nrow(anova) == 1) return(c(MSE = anova$`Mean Sq`[nrow(anova)]))
  MSE = anova$`Mean Sq`[nrow(anova)]
  table = data.frame(anova[,c('F value','Pr(>F)')])
  table = table[-nrow(table),,drop=FALSE]
  names(table) = c('Fvalue','Pvalue')
  c(MSE=MSE,unlist(lapply(names(table),function(x) {value = table[[x]];names(value) = paste(rownames(table),x,sep='::');value})))
}


#' Title
#'
#' @param Y n x p matrix with each column to be tested separately against X_matrices
#' @param X_matrices n x bm matrix of m nxb design matrices each with b columns
#'   to be jointly tested
#' @param X_cov n x a matrix of covariates constant for each test
#' @param cV matrix R such that R'R = V. Cholesky decomposition of the nxn
#'   covariance matrix V
#' @param testGroups vector of length bm with a uniq value for each of the m
#'   tests repeated b times.
#' @param assign optional vector of length b that breaks up the b column of each
#'   design matrix into multiple groups for the ANOVA
#' @param ML currently not used
#' @param REML currently not used
#' @param mc.cores number of cores to use for parallel solutions
#' @param verbose should progress be reported?
#'
#' @return list with elements matrices of test statistics for each test
#' @export
#'
#' @examples
EMMAX_ANOVA = function(formula,data,markers,genotypeID,cholL_Sigma_inv,mc.cores = RcppParallel::defaultNumThreads()-1,verbose = T) {
  require(foreach)
  require(doParallel)
  # step 1: parse formula, extract the terms with "X" meaning marker, create two formulas, one for X_cov the other for X_base
  # X_cov: covariates that are constant for each test
  # X_base: model matrix for marker effects: diag(marker[,j]) %*% X_base is the design matrix for marker j

  terms_full = terms(formula)
  factors = attr(terms_full,'factors')
  if(!'X' %in% rownames(factors) ) stop('no marker terms "X" found in formula')
  if(length(grep('X',colnames(factors[,factors['X',]==0,drop=FALSE])))) stop('character "X" is a variable name. Please re-name')

  X_terms = factors['X',]>0
  if(length(X_terms) == 0) stop('no marker terms "X" found in formula')
  base_formula = update(formula,sprintf('~.-%s',paste(attr(terms_full,'term.labels')[X_terms],collapse='-')))

  marker_formula = update(formula,sprintf('~.-%s-1',paste(c('1',attr(terms_full,'term.labels')[!X_terms]),collapse='-')))

  mf = model.frame(base_formula,data)
  if(!attr(attr(mf,'terms'),'response')) stop('no response specified')
  Y = mf[[1]]
  X_cov = model.matrix(base_formula,mf)
  qr_X_cov = qr(X_cov)
  X_cov = X_cov[,qr_X_cov$pivot[1:qr_X_cov$rank]]


  # given a data matrix Y, matrix of covariates (constant across tests), a set
  # of X matrices concatenated together, and cV = chol(V): calculates b_hat,
  # SE(b_hat), and the F-test that all b_hat == 0 for the tests for each column
  # of Y, projecting out X_cov


  if(is(Y,'vector')) Y = matrix(Y)
  if(is.null(colnames(Y))) colnames(Y) = 1:ncol(Y)
  # recover()
  # P matrix to project out X_cov
  n = nrow(Y)
  b_cov = ncol(X_cov)

  cVi_Y = as.matrix(cholL_Sigma_inv %**% Y)
  cVi_Xcov = as.matrix(cholL_Sigma_inv %**% X_cov)
  t_cVi_Xcov = t(cVi_Xcov)
  Xcovt_Vi_Xcov = crossprod(cVi_Xcov)
  Xcovt_Vi_Xcov_inv = solve(Xcovt_Vi_Xcov)
  cVi_Xcov_Xcovt_Vi_Xcov_inv = cVi_Xcov %**% Xcovt_Vi_Xcov_inv

  PY = cVi_Y - cVi_Xcov_Xcovt_Vi_Xcov_inv %**% (t_cVi_Xcov %**% cVi_Y)


  if(verbose) print('Done setup')
  registerDoParallel(mc.cores)
  chunks = unique(c(seq(0,ncol(markers),by = 1000*mc.cores),ncol(markers)))
  results = foreach(j = 2:length(chunks),.combine = c) %do% {
    index = seq(chunks[j-1]+1,chunks[j])
    if(verbose) print(sprintf('%d of %d',index[1],ncol(markers)))
    foreach(i = index) %dopar% {
      data$X = markers[data[[genotypeID]],i]
      data$X = data$X - as.numeric(names(sort(table(data$X),decreasing = T)[1]))
      X_design = Matrix::sparse.model.matrix(marker_formula,data)
      assign = attr(X_design,'assign')
      X_design = Matrix::drop0(X_design)
      # find non-zero rows of X. Generally should be a large percentage if most minor alleles are rare
      j = Matrix::rowSums(X_design != 0) != 0

      # rotate X_design:
      cVi_design = partial_matrix_multiply_toDense(cholL_Sigma_inv,X_design,j)
      cVi_design_sparse = as(cVi_design,'dgCMatrix')
      PX_design = cVi_design - cVi_Xcov_Xcovt_Vi_Xcov_inv %*% (t_cVi_Xcov %**% cVi_design_sparse)
      colnames(PX_design) = colnames(X_design)

      object = lm.fit(PX_design,PY)

      object$df.residual = object$df.residual - b_cov
      object$terms = terms(marker_formula)
      object$assign = assign
      class(object) = 'lm'
      p1 = 1L:object$rank
      if(object$rank>0) {
        sqrtR <- sqrt(diag(chol2inv(object$qr$qr[p1, p1, drop = FALSE])))
      } else{
        sqrtR = NA
      }
      # Need to repeat these lines separately for each trait if Y is a matrix, otherwise we don't get the univariate anova results
      anovas = do.call(rbind,lapply(1:ncol(PY),function(t) {
        lm_object = list(
          residuals = as.matrix(object$residuals)[,t],
          fitted.values = as.matrix(object$fitted.values)[,t],
          effects = as.matrix(object$effects)[,t],
          df.residual = object$df.residual,
          rank = object$rank,
          terms = object$terms,
          assign = object$assign,
          qr = object$qr
        )
        class(lm_object) = 'lm'
        anova = anova(lm_object)
        # print(anova)
        # print(anova_table(anova))
        anova_table(anova)
      }))
      rownames(anovas) = colnames(Y)
      MSEs = anovas[,'MSE']
      beta_hats = as.matrix(object$coefficients)[p1,,drop=FALSE]
      SEs = outer(sqrtR,MSEs)
      rownames(SEs) = rownames(as.matrix(object$coefficients))[p1]
      colnames(beta_hats) = colnames(SEs) = colnames(Y)
      results = list(
        anova = anovas,
        beta_hats = beta_hats,
        SEs = SEs
      )
    }
  }
  collect_results(results)
}



