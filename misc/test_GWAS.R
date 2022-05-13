library(JointGWAS)
d = data.frame(ID = 1:100)
d$X = sample(c(0,1),nrow(d),replace=T,prob = c(.9,.1))
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

m = relmatLmer(y~0+Trait+X:Trait+(0+Trait|ID),d_tall,relmat = list(ID = K))
vars = as.data.frame(VarCorr(m))$vcov

Ghat = matrix(c(vars[1],vars[3],vars[3],vars[2]),2)
Rhat = diag(vars[4],2)
colnames(Ghat) = colnames(Rhat) = rownames(Ghat) = rownames(Rhat) = c(1,2)

cholL_Sigma_inv = make_cholL_Sigma_inv(d_tall,'y','ID','Trait',list(list(Row=K,Column=Ghat),list(Column=Rhat)))

markers = matrix(d$X,nrow = nrow(d),ncol=1:2);rownames(markers) = d$ID
res = fastEMMAX_ANOVAs(formula=y~0+Trait+X:Trait,d_tall,markers,'ID',cholL_Sigma_inv,1)
res
summary(m)$coef[3:4,1] - res[[1]]$beta_hats[1,]
anova(m)$F[-1] - res[[1]]$anova[1,2]
