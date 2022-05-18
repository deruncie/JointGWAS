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
  table = data.frame(anova[,c('Df','F value','Pr(>F)')])
  table = table[-nrow(table),,drop=FALSE]
  names(table) = c('Df','Fvalue','Pvalue')
  c(MSE=MSE,unlist(lapply(names(table),function(x) {value = table[[x]];names(value) = paste(rownames(table),x,sep='::');value})))
}


#' Joint ANOVA GWAS for multiple correlated traits
#'
#' Runs a Genome Wide Association study for multiple correlated traits, performing an ANOVA for
#' each marker on the joint effect across all traits.
#'
#' @param formula lm-style formula for the model to be fit for each marker.
#' The formula must include the special variable "X" which represents a marker.
#' \code{data} MUST NOT include "X" as a column, as this will be created in \code{data}
#' for each marker in the GWAS, replacing anything in that slot. Typically the formula
#' will include a main effect of each trait and an interaction between the marker and the trait.
#' @param data data.frame with at a minimum columns for trait and \code{genotypeID} to be used to match
#' up with \code{markers}. MUST NOT include a column titled "X".
#' @param markers matrix of marker genotypes with individuals as rows and rownames corresponding to the
#' \code{genotypeID} column of \code{data}
#' @param genotypeID column of \code{data} representing the genotype identifiers (rownames of \code{markers})
#' @param cholL_Sigma_inv lower-triangular matrix of the inverse of the Cholesky decomposition of the full variance-covariance
#' matrix of the data. Typically calculated with \code{make_cholL_Sigma_inv}
#' @param mc.cores Number of processor cores to use for parallel calculations
#' @param verbose Should progress across markers be reported?
#' @param MAC_filter minimum minor allele counts per trait
#' @param MAF_filter minimum minor allele frequency per trait
#'
#' @return list of results with three elements:
#'     anovas: MSE, F-values and p-values for each term involving "X" for each marker
#'     beta_hats: coefficient estimates for each X term
#'     SEs: standard errors for each X term coefficient.
#' @export
#'
EMMAX_ANOVA = function(formula,data,markers,genotypeID,cholL_Sigma_inv,mc.cores = RcppParallel::defaultNumThreads()-1,verbose = T, MAC_filter = NULL, MAF_filter = NULL) {
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


  data$X = 1
  X_design_base = Matrix::sparse.model.matrix(marker_formula,data)!=0
  n_per_coef = Matrix::colSums(X_design_base != 0)

  if(verbose) print('Done setup')
  registerDoParallel(mc.cores)
  chunks = unique(c(seq(0,ncol(markers),by = 1000*mc.cores),ncol(markers)))
  results = foreach(j = 2:length(chunks),.combine = c) %do% {
    index = seq(chunks[j-1]+1,chunks[j])
    if(verbose) print(sprintf('%d of %d',index[1],ncol(markers)))
    foreach(i = index) %dopar% {
      # recover()
      data$X = markers[data[[genotypeID]],i]
      data$X = data$X - as.numeric(names(sort(table(data$X),decreasing = T)[1]))
      X_design = Matrix::sparse.model.matrix(marker_formula,data)
      assign = attr(X_design,'assign')
      X_design = Matrix::drop0(X_design)
      if(!is.null(MAF_filter) || !is.null(MAC_filter)) {
        macs = Matrix::colSums(X_design != 0)
        drop_cols = rep(F,ncol(X_design))
        if(!is.null(MAC_filter)) drop_cols[macs < MAC_filter] = T
        if(!is.null(MAF_filter)) drop_cols[macs/n_per_coef < MAF_filter] = T
        if(any(drop_cols)) {
          X_design = X_design[,!drop_cols,drop=FALSE]
          assign = assign[!drop_cols,drop=FALSE]
        }
      }
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
          effects = if(!is.null(object$effects)) as.matrix(object$effects)[,t],
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
      if(length(object$coefficients)>0) {
        beta_hats = as.matrix(object$coefficients)[p1,,drop=FALSE]
      } else {
        beta_hats = matrix(NA,ncol = ncol(Y))
      }
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
  output = collect_results(results)
  if(length(output) == 1) output = output[[1]]
  return(output)
}



