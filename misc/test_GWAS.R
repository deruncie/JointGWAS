library(JointGWAS)
n=100
d = data.frame(ID = 1:n)
d$X = sample(c(0,1),nrow(d),replace=T,prob = c(.9,.1))
d$env = rnorm(nrow(d))
K = tcrossprod(MegaLMM::rstdnorm_mat(nrow(d),nrow(d))) + diag(1,nrow(d))
K = K/mean(diag(K))
rownames(K) = colnames(K) = d$ID

G = .5+diag(.2,2)

U = t(chol(K)) %*% matrix(rnorm(nrow(d)*2),ncol=2) %*% chol(G)
E = matrix(rnorm(nrow(d)*2),ncol=2)
Y = U=E

d_tall = rbind(data.frame(d,y=Y[,1],Trait=1),data.frame(d,y=Y[,2],Trait=2))
d_tall$Trait = factor(d_tall$Trait)

library(lme4qtl)
library(lme4)

m = relmatLmer(y~0+Trait+Trait:env + X:Trait+X:Trait:env + (0+Trait|ID),d_tall,relmat = list(ID = K))
vars = as.data.frame(VarCorr(m))$vcov

Ghat = matrix(c(vars[1],vars[3],vars[3],vars[2]),2)
Rhat = diag(vars[4],2)
colnames(Ghat) = colnames(Rhat) = rownames(Ghat) = rownames(Rhat) = c(1,2)

markers = matrix(d$X,nrow = nrow(d),ncol=1:2);rownames(markers) = d$ID
cholL_Sigma_inv = make_cholL_Sigma_inv(d_tall,'y','ID','Trait',list(list(Row=K,Column=Ghat),list(Column=Rhat)))
res = EMMAX_ANOVA(formula=y~0+Trait+Trait:env + X:Trait+X:Trait:env,d_tall,markers,'ID',cholL_Sigma_inv,1)
sK = svd(K)
sK2 = simultaneous_diagonalize(K,diag(1,n))
sGR = simultaneous_diagonalize(Ghat,Rhat)
res2 = EMMAX_ANOVA_matrix(Y~X+env+X:env,d,markers,'ID',svd_matrices = list(sK,sGR),1)
res
res2
summary(m)$coef
anova(m)

d_tall2 = d_tall[-c(1:(n/2)),]
cholL_Sigma_inv = make_cholL_Sigma_inv(d_tall2,'y','ID','Trait',list(list(Row=K,Column=Ghat),list(Column=Rhat)))
markers = do.call(cbind,lapply(1:10000,function(i) rbinom(n,1,runif(1,0,.25))));rownames(markers) = d$ID
markers = markers[,colSums(markers)>=10 & colSums(1-markers) >= 10]
res = EMMAX_ANOVA(formula=y~0+Trait+Trait:env + X:Trait+X:Trait:env,d_tall2,markers,'ID',cholL_Sigma_inv,1,MAF_filter = 0.05)
qq(subset(res$anova,Trait.X..Df == 1)[,6])
qq(subset(res$anova,Trait.X..Df == 2)[,6])
res2 = EMMAX_ANOVA(formula=y~0+Trait+Trait:env + X:Trait+X:Trait:env,d_tall2,markers,'ID',cholL_Sigma_inv,1,MAC_filter = 0)
qq(res2$anova$Trait.X..Pvalue[res$anova$Trait.X..Df==1])
qq(res2$anova$Trait.X..Pvalue[res$anova$Trait.X..Df==2])
# testing with more traits and cis_markers
t = 10
n = 100
Y = MegaLMM::rstdnorm_mat(n,t)
d = data.frame(ID = 1:n)
d$X = sample(c(0,1),nrow(d),replace=T,prob = c(.9,.1))
d$env = rnorm(nrow(d))
K = tcrossprod(MegaLMM::rstdnorm_mat(nrow(d),nrow(d))) + diag(1,nrow(d))
K = K/mean(diag(K))
K = diag(1,n)
rownames(K) = colnames(K) = d$ID
G = tcrossprod(MegaLMM::rstdnorm_mat(t,t)) + diag(1,t)
R = tcrossprod(MegaLMM::rstdnorm_mat(t,t)) + diag(1,t)
G = R = diag(.5,t)
sK = svd(K)
sGR = simultaneous_diagonalize(G,R)
markers = matrix(d$X,nrow = nrow(d),ncol=1:2);rownames(markers) = d$ID
cis_markers = list()
cis_markers[[1]] = c(5,7,9)
res2 = EMMAX_ANOVA_matrix(Y~X,d,markers,'ID',svd_matrices = list(sK,sGR),mc.cores = 1,cis_markers=cis_markers)
res2
d_tall = data.frame(Trait = factor(rep(1:t,each=n)),X = d$X[rep(1:n,t)],y = c(Y))
anova(lm(y~Trait+X:Trait,d_tall))
X1 = model.matrix(~Trait,d_tall)
X2 = model.matrix(~X:Trait-1,d_tall)
anova(lm(y~1+X1+X2[,cis_markers[[1]]]+X2[,-cis_markers[[1]]],d_tall))

#
# t(sK$u) %*% K %*% sK$u
# t(sGR$S) %*% Ghat %*% sGR$S
# t(sGR$S) %*% Rhat %*% sGR$S
# sGR$d
# X = model.matrix(~X,d)
#
# D = 1 + c(outer(sK$d,sGR$d))
# rY = 1/sqrt(D)*(t(sK$u) %*% Y %*% sGR$S)
# sX = t(sK$u) %*% X
# rX = 1/sqrt(D)*(do.call(cbind,lapply(1:ncol(sX),function(i) matrix(outer(sX[,i],sGR$S),ncol=2))))
# rX2 = 1/sqrt(D)*matrix(outer(sX,sGR$S),ncol=2)
# # rX = (do.call(cbind,lapply(1:ncol(sX),function(i) c(sX[,rep(i,2)] %*% sGR$S))))
# # rX = (do.call(cbind,lapply(1:ncol(sX),function(i) matrix(outer(sX[,i],sGR$S),ncol=2))))
#
# # NOTE: solve(crossprod(rX)) is not block-diagonal by trait, so will have to fit all traits jointly
# # NOTE: but above code for forming rX should be useful to avoid the tn x tn operation? Unless the outer call is tn x tn (I think it's tn x t)
#
# summary(lm(c(rY)~0+rX))
# summary(m)$coef
#
# s = kronecker(t(sGR$S),t(sK$u))
# y = c(Y)
# ry = 1/sqrt(D)*(s %*% y)
# plot(ry,c(rY))
# X2 = cbind(bdiag(X[,1],X[,1]),bdiag(X[,2],X[,2]))
# X2 = bdiag(X,X)
# # rX2 = 1/sqrt(D)*(s %*% X2)
# rX2 = (s %*% X2)
# plot(rX2[,1],rX[,1])
#
# Sigma = kronecker(Ghat,K) + kronecker(Rhat,diag(1,nrow(Y)))
# diag(s %*% Sigma %*% t(s) - diag(kronecker(diag(sGR$d),diag(sK$d))))
# diag(kronecker(diag(sGR$d),diag(sK$d))) - c(outer(sK$d,sGR$d))
# diag(s %*% Sigma %*% t(s)) - 1 - c(outer(sK$d,sGR$d))
#
# f = function(formula,data,svd_matrices = list(sK,sGR)) {
#   recover()
# }
# f(Y~X,d)

