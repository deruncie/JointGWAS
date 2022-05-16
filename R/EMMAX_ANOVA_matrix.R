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

EMMAX_ANOVA_matrix = function(formula,data,markers,genotypeID,svd_matrices,mc.cores = RcppParallel::defaultNumThreads()-1,verbose = T) {
  # svd_matrices is a list with two components:
  #   svd_GR = simultaneous_diagonalize(R,Ginvsq)
  #   svd_K = svd(K)
  #   We have Sigma = kronecker(K,G) + kronecker(I,R)
  #   If we pre-multiply by kronecker(t(svd_K$u),t(svd_GR$S)), the result is a diagonal matrix
  #   We can de-correlate Y by t(svd_K$u) %*% Y %*% svd_GR$S
  #   We can de-correlate rows of X by t(svd_K$u) %*% X
  #   How to de-correlated columns of B? B %*% svd_GR$S  How to then test the columns of B?

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
  X_cov = X_cov[,qr_X_cov$pivot[1:qr_X_cov$rank],drop=FALSE]

  if(is(Y,'vector')) Y = matrix(Y)
  if(is.null(colnames(Y))) colnames(Y) = 1:ncol(Y)
  n = nrow(Y)
  t = ncol(Y)
  b_cov = ncol(X_cov)


  # Extract rotation matrices
  rotate_rows = t(svd_matrices[[1]]$u)
  rotate_cols = svd_matrices[[2]]$u
  diag_vars = 1 + c(outer(svd_matrices[[1]]$d,svd_matrices[[2]]$d))

  cVi_Y = as.matrix(1/sqrt(diag_vars)*c(rotate_rows %*% Y %*% rotate_cols))

  rX_cov = rotate_rows %*% X_cov
  cVi_Xcov = 1/sqrt(diag_vars) * do.call(cbind,lapply(1:ncol(X_cov),function(i) matrix(outer(rX_cov[,i],rotate_cols),ncol=t)))

  t_cVi_Xcov= t(cVi_Xcov)
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
      rX_design = partial_matrix_multiply_toDense(rotate_rows,X_design,j)
      cVi_design = 1/sqrt(diag_vars) * do.call(cbind,lapply(1:ncol(rX_design),function(i) matrix(outer(rX_design[,i],rotate_cols),ncol=t)))
      PX_design = cVi_design - cVi_Xcov_Xcovt_Vi_Xcov_inv %*% (t_cVi_Xcov %**% cVi_design)
      colnames(PX_design) = paste(colnames(Y)[rep(1:t,each=ncol(X_design))],colnames(X_design)[rep(1:ncol(X_design),t)],sep='::')

      object = lm.fit(PX_design,PY)

      object$df.residual = object$df.residual - b_cov*t
      object$terms = terms(marker_formula)
      object$assign = rep(assign,each = t)
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
      rownames(anovas) = '1';
      MSEs = anovas[,'MSE']
      beta_hats = as.matrix(object$coefficients)[p1,,drop=FALSE]
      SEs = outer(sqrtR,MSEs)
      rownames(SEs) = rownames(as.matrix(object$coefficients))[p1]
      colnames(beta_hats) = colnames(SEs) = '1'#colnames(Y)
      results = list(
        anova = anovas,
        beta_hats = beta_hats,
        SEs = SEs
      )
    }
  }
  collect_results(results)[[1]]
}
