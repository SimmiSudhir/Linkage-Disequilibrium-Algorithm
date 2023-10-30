## cross validation using continuous model
## usign asreml using continuous genotypes 58K

rm(list=ls())

library(data.table)
library(asreml)

Markers <- as.data.frame(fread("M_con58_1317.txt",header=TRUE))
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
fwrite(Ginv.sparse,file='GRM_1317_con58.txt',
       row.names=F, col.names=T,sep= " ")

########################################################
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
#Pheno1 <- cbind(clones, Pheno$CCS)
Pheno1 <- cbind(clones, Pheno$Fibre)

colnames (Pheno1)[1] <- "Clone"
colnames(Pheno1)[2] <- "TCH"
head(Pheno1)
str(Pheno1)
Pheno1$Clone <- as.factor(Pheno1$Clone)

ahatinv<-read.table("GRM_1317_con58.txt",h=T)

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
# > apply(Tab[-1],2, function(x)c(mean(x),sd(x)))
# RMSE        ACC
# [1,] 8.7670998 0.41916716
# [2,] 0.3405316 0.06243994
# > Tab
# PT     RMSE       ACC
# 1  1 8.792405 0.4420247
# 2  2 8.944877 0.3658962
# 3  3 8.342784 0.5119747
# 4  4 9.214678 0.3590308
# 5  5 8.540756 0.4169093
# > h2_cv
# [,1]
# [1,] 0.3768402
# [2,] 0.4193458
# [3,] 0.3390334
# [4,] 0.3662916
# [5,] 0.3766081

# ## For CCS
# > apply(Tab[-1],2, function(x)c(mean(x),sd(x)))
# RMSE        ACC
# [1,] 0.56279439 0.46116404
# [2,] 0.04365198 0.03924793
# > Tab
# PT      RMSE       ACC
# 1  1 0.5850029 0.4244641
# 2  2 0.4918600 0.4961522
# 3  3 0.6045536 0.4167417
# 4  4 0.5534524 0.5008597
# 5  5 0.5791031 0.4676026
# > h2_cv
# [,1]
# [1,] 0.5176348
# [2,] 0.4985151
# [3,] 0.5368264
# [4,] 0.5470563
# [5,] 0.5146957

# ## For Fibre
# apply(Tab[-1],2, function(x)c(mean(x),sd(x)))
# RMSE        ACC
# [1,] 1.07224661 0.48685828
# [2,] 0.05066751 0.03888538
# > Tab
# PT     RMSE       ACC
# 1  1 1.095290 0.4871912
# 2  2 1.080214 0.5139845
# 3  3 1.034373 0.4970273
# 4  4 1.011583 0.4206048
# 5  5 1.139773 0.5154837
# > h2_cv
# [,1]
# [1,] 0.5201569
# [2,] 0.5462467
# [3,] 0.5357063
# [4,] 0.5930602
# [5,] 0.5410663
