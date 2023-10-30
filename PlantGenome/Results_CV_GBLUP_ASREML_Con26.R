## cross validation using continuous model
## usign asreml using continuous genotypes 26K

rm(list=ls())

library(data.table)
library(asreml)

Markers <- as.data.frame(fread("M_con_1317.txt",header=TRUE))
dim(Markers) ## 1317 25753
Markers[1:5,1:5]

Clones_1317 <- as.data.frame(fread("1317_clones.txt", header=TRUE))
head(Clones_1317)
colnames(Clones_1317)[1] <- "Clone"
dim(Clones_1317) ## 1317 by 1

####

rownames(Markers) <- Clones_1317$Clone
Markers[1:5,1:5]

X <- scale(Markers, center=TRUE, scale=TRUE)
m= nrow(X)
dim(X) ## 1317 25753
Z1 <- (X %*% t(X))
dim(Z1) ## 1317 1317

G = Z1/ncol(X)
G[1:5,1:5]

####Check the diagonal sum

mean(diag(G)) ## 0.9992407
det(G) ## 0
##### Check the mean of off-diagonal
mean(G[row(G)!=col(G)]) ## -0.0007593014

G <- G+diag(m)*0.001
Ginv <- solve(G) ## 1318 1318
Ginv[1:5,1:5]
dim(Ginv)

###### inverse of additive relationship matrix, saved in sparse form ##########

source('C:/Users/uqsyada1/OneDrive - The University of Queensland/PhD/Trainings and workshop/201910_AsREML/Asreml-R workshop/Asreml(4-5)/Material_UQ45/full2sparse.R')
Ginv.sparse<-full2sparse(Ginv)
colnames(Ginv.sparse)<-c('Row','Column','Ainverse')
head(Ginv.sparse)
dim(Ginv.sparse)
fwrite(Ginv.sparse,file='GRM_1317_con.txt',
       row.names=F, col.names=T,sep= " ")


#################################################
## Continuous Matrix with 26K markers 

library(asreml)
library(data.table)
rm(list=ls())
# Reading Phenotypic data
Pheno <- as.data.frame(fread("Blues_1317.txt", header=T))
head(Pheno)
clones <- read.table ("1317_clones.txt", header = T)
colnames(clones)[1] <- "x"


#Pheno1 <- cbind(clones, Pheno$TCH)
# Pheno1 <- cbind(clones, Pheno$CCS)
Pheno1 <- cbind(clones, Pheno$Fibre)

colnames (Pheno1)[1] <- "Clone"
colnames(Pheno1)[2] <- "TCH"
head(Pheno1)
str(Pheno1)
Pheno1$Clone <- as.factor(Pheno1$Clone)

ahatinv<-read.table("GRM_1317_con.txt",h=T)

ahatinv<-as.matrix.data.frame(ahatinv)
attr(ahatinv,"rowNames")<-as.character(clones$x)
attr(ahatinv,"colNames")<-as.character(clones$x)
attr(ahatinv,"INVERSE")<-TRUE
head(ahatinv)
library(asreml)


n <- length(Pheno1$Clone)
K <- 5
set.seed(123)
group <- sample(c(1:K),n, replace=TRUE)
table(group)
h2_cv <- matrix(data=NA, nrow=K, ncol=1)
Tab <- data.frame(PT=1:K, RMSE=NA, ACC=NA)

for(g in 1:K){
  test_idx <- group == g
  train_idx <- !test_idx
  y_train <- Pheno1$TCH[train_idx]
  y_test <- Pheno1$TCH[test_idx]
  y <- Pheno1$TCH
  ycv <- y
  for(j in 1:n) {
    if (group[j]==g) {ycv[j] <- NA}
  }
  modelGBLUP_cv = asreml(fixed=ycv ~1,
                         random=~vm(Clone,ahatinv),
                         #random=~Clone,
                         workspace=128e06,na.action = na.method(y="include"),
                         data= Pheno1)
  
  predGBLUPcv <- predict(modelGBLUP_cv, classify="Clone", sed=T)$pval
  colnames(predGBLUPcv)[2] <- "Solution"
  head(predGBLUPcv)
  Yp_ts <- predGBLUPcv[test_idx, 2]
  
  h2_GBLUPcv <- vpredict(modelGBLUP_cv, h2~((V1/(V1+V2))))
  h2_cv[g] <- as.numeric(h2_GBLUPcv[1])
  
  Tab$RMSE[g] <- sqrt(mean((y_test - Yp_ts)^2, na.rm=TRUE))
  Tab$ACC[g] <- cor(y_test, Yp_ts, method = "pearson", use = "complete.obs")
}
apply(Tab[-1],2, function(x)c(mean(x),sd(x)))

# ## For TCH 
# RMSE        ACC
# [1,] 8.7958442 0.41253949
# [2,] 0.3388486 0.06377004
# > Tab
# PT     RMSE       ACC
# 1  1 8.838840 0.4322728
# 2  2 8.960676 0.3581107
# 3  3 8.374853 0.5087609
# 4  4 9.241766 0.3520839
# 5  5 8.563088 0.4114692
# > h2_cv
# [,1]
# [1,] 0.3467495
# [2,] 0.3911750
# [3,] 0.3055542
# [4,] 0.3448930
# [5,] 0.3562646
# ## For CCS
# > apply(Tab[-1],2, function(x)c(mean(x),sd(x)))
# RMSE       ACC
# [1,] 0.56639774 0.4501782
# [2,] 0.04485325 0.0425071
# > Tab
# PT      RMSE       ACC
# 1  1 0.5906454 0.4062577
# 2  2 0.4931203 0.4925634
# 3  3 0.6065214 0.4104513
# 4  4 0.5560157 0.4939652
# 5  5 0.5856858 0.4476535
# > h2_cv
# [,1]
# [1,] 0.4984861
# [2,] 0.4680380
# [3,] 0.5128553
# [4,] 0.5233561
# [5,] 0.4901299

# ## For Fibre
# > apply(Tab[-1],2, function(x)c(mean(x),sd(x)))
# RMSE        ACC
# [1,] 1.08327905 0.47039962
# [2,] 0.05523692 0.03725541
# > Tab
# PT     RMSE       ACC
# 1  1 1.104319 0.4728117
# 2  2 1.089488 0.5014148
# 3  3 1.042877 0.4850445
# 4  4 1.019130 0.4062751
# 5  5 1.160582 0.4864520
# > h2_cv
# [,1]
# [1,] 0.4923954
# [2,] 0.5125556
# [3,] 0.5055057
# [4,] 0.5661173
# [5,] 0.5137050
# > 

