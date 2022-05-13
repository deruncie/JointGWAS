collect_results = function(results_list) {
  results = list(
    beta_hats = do.call(rbind,lapply(results_list,function(x) x$beta_hats)),
    SEs = do.call(rbind,lapply(results_list,function(x) x$SE)),
    RSSs = do.call(rbind,lapply(results_list,function(x) x$RSS)),
    sigmas = do.call(rbind,lapply(results_list,function(x) x$sigma)),
    Fvalues = do.call(rbind,lapply(results_list,function(x) x$Fvalue)),
    Pvalues = do.call(rbind,lapply(results_list,function(x) x$p_value)),
    V_log_dets = do.call(c,lapply(results_list,function(x) x$V_log_det)),
    V_star_inv_log_dets = do.call(c,lapply(results_list,function(x) x$V_star_inv_log_det)),
    MLs = do.call(rbind,lapply(results_list,function(x) x$ML)),
    REMLs = do.call(rbind,lapply(results_list,function(x) x$REML))
  )
  results
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
fastEMMAX_F_tests = function(Y,X_matrices,X_cov,cV,testGroups,assign = NULL,ML=FALSE,REML=FALSE,mc.cores = RcppParallel::defaultNumThreads()-1,verbose = F) {
  # given a data matrix Y, matrix of covariates (constant across tests), a set
  # of X matrices concatenated together, and cV = chol(V): calculates b_hat,
  # SE(b_hat), and the F-test that all b_hat == 0 for the tests for each column
  # of Y, projecting out X_cov

  require(foreach)
  require(doParallel)
  require(MegaLMM)
  determinants = ML || REML
  if(is(Y,'vector')) Y = matrix(Y)
  if(is.null(assign)) assign = rep(1,table(testGroups)[1])
  # recover()
  # P matrix to project out X_cov
  n = nrow(Y)
  b_cov = ncol(X_cov)

  # if(verbose) print('Y')
  cVi_Y = forwardsolve(t(cV),Y)
  # if(verbose) print('X_cov')
  cVi_Xcov = forwardsolve(t(cV),X_cov)
  Xcovt_Vi_Xcov = crossprod(cVi_Xcov)
  sXcovt_Vi_Xcov = svd(Xcovt_Vi_Xcov)
  tol = 1e-7
  Positive = sXcovt_Vi_Xcov$d > tol
  b_cov = sum(Positive)
  Xcovt_Vi_Xcov_inv = sXcovt_Vi_Xcov$v[,Positive,drop=FALSE] %*% ((1/sXcovt_Vi_Xcov$d[Positive])*t(sXcovt_Vi_Xcov$u[,Positive,drop=FALSE]))
  PY = cVi_Y - cVi_Xcov%**% (Xcovt_Vi_Xcov_inv %**% (t(cVi_Xcov) %**% cVi_Y))

  if(REML) {
    log_det_Xcovt_Vi_Xcov = determinant(Xcovt_Vi_Xcov,logarithm = TRUE)$modulus
    inv_Xcovt_Xcov = solve(crossprod(X_cov))
    log_det_Xcovt_Xcov = -determinant(inv_Xcovt_Xcov,logarithm = TRUE)$modulus
  }

  # joint rotations of all columns of X_matrices
  if(verbose) print('X_matrices')
  cVi_X_matrices = forwardsolve(t(cV),X_matrices)
  PXni = cVi_X_matrices - cVi_Xcov %**% (Xcovt_Vi_Xcov_inv %**% (t(cVi_Xcov) %**% cVi_X_matrices))

  rm('cVi_X_matrices')
  gc()

  if(verbose) print('Done setup')
  # registerDoParallel(mc.cores)
  results = foreach(test=split(1:length(testGroups),testGroups)) %do% {
    b_hat = rep(NA,length(test))
    SE = rep(NA,length(test))
    SSE = NA
    s2_hat = NA
    Fvalue = NA
    p_value = NA

    X_test_i = PXni[,test,drop=FALSE]
    test_sub = colSums(X_test_i^2) > 1e-10
    asgn = assign[test_sub]
    if(sum(test_sub)>0) {
      object = lm.fit(X_test_i[,test_sub,drop=FALSE],PY)
      ssr = sum(object$residuals^2)
      mss = sum(object$fitted.values^2)
      dfr <- df.residual(object) - b_cov
      p <- object$rank
      p1 <- 1L:p
      comp <- object$effects[p1]
      ss <- c(vapply(split(comp^2, asgn), sum, 1), ssr)
      df <- c(lengths(split(asgn, asgn)), dfr)
      ms <- (ss/df)[max(asgn)]
      Fvalue <- ms/(ssr/dfr)
      p_value <- pf(Fvalue, df[max(asgn)], dfr, lower.tail = FALSE)

      s2_hat = ssr/dfr
      R <- chol2inv(object$qr$qr[p1, p1, drop = FALSE])
      b_hat[test_sub] = object$coefficients
      SE[test_sub] = sqrt(diag(R) * s2_hat)
    }
    result = list(beta_hats = b_hat,
                  SE = SE,
                  RSS = SSE,
                  sigma = sqrt(s2_hat),
                  Fvalue = Fvalue,
                  p_value = p_value
    )
    result
  }
  collect_results(results)
}
