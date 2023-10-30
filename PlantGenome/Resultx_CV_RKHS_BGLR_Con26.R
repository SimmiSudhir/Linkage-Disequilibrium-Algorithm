############ Multi Kernel RKHS regression
## using continuous 26K marker
rm(list=ls())

library(data.table)
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

## genotypes
geno <- as.data.frame(fread("M_con_1317.txt",header=TRUE))
dim(geno) ## 1317 by 25753
geno[1:5,1:5]

clones_1317 <- as.data.frame(fread("1317_clones.txt",header=TRUE))
head(clones_1317)
colnames(clones_1317)[1] <- "Clone"

rownames(geno) <- clones_1317$Clone
geno[1:5,1:5]

## scale the marker matrix
library(BGLR)
X <- scale(geno, center=TRUE, scale=TRUE)
dim(X) ## 1317 25753
p <- ncol(X)

## calculate euclidean matrix #######

D <- (as.matrix(dist(X,method='euclidean'))^2)/p
D[1:5,1:5]
dim(D)
######## write the distance matrix in directory

fwrite(D,file="DM_1317_con.txt",row.names=F,col.names=T, sep=" ")
## Read the file from the directory ############################

dis_mat <- as.data.frame(fread("DM_1317_con.txt",header=TRUE))
dis_mat <- as.matrix(dis_mat)
dim(dis_mat)
dis_mat[1:5,1:5]

h<-0.5*c(1/5,1,5)
#3# Kernel Averaging using BGLR
ETA<-list(list(K=exp(-h[1]*dis_mat),model='RKHS'),
          list(K=exp(-h[2]*dis_mat),model='RKHS'),
          list(K=exp(-h[3]*dis_mat),model='RKHS'))

## 5 fold cross validation 

n <- length(Pheno1$Clone)
K <- 5
set.seed(123)
group <- sample(c(1:K),n, replace=TRUE)
table(group)
#h2_cv <- matrix(data=NA, nrow=K, ncol=1)
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
  modelBGLR_cv = BGLR(y=ycv, ETA=ETA,nIter=5000, burnIn=2000)
  
  
  Yp_ts <- modelBGLR_cv$yHat[test_idx]
  
  #h2_GBLUPcv <- vpredict(modelGBLUP_cv, h2~((V1/(V1+V2))))
  #h2_cv[g] <- as.numeric(h2_GBLUPcv[1])
  
  Tab$RMSE[g] <- sqrt(mean((y_test - Yp_ts)^2, na.rm=TRUE))
  Tab$ACC[g] <- cor(y_test, Yp_ts, method = "pearson", use = "complete.obs")
}
apply(Tab[-1],2, function(x)c(mean(x),sd(x)))

## 
## For TCH
#         RMSE        ACC
# [1,] 8.7080351 0.43053730
# [2,] 0.6203632 0.05891421
# > Tab
# PT     RMSE       ACC
# 1  1 8.119215 0.4430472
# 2  2 8.387009 0.4705296
# 3  3 9.290129 0.3865206
# 4  4 8.281534 0.4975483
# 5  5 9.462289 0.3550409

## For CCS 
#         RMSE        ACC
# [1,] 0.56632076 0.45919880
# [2,] 0.02451068 0.05628345
# > Tab
# PT      RMSE       ACC
# 1  1 0.5664826 0.3948976
# 2  2 0.5284454 0.4906634
# 3  3 0.5638487 0.4881188
# 4  4 0.5776111 0.5190475
# 5  5 0.5952161 0.4032667

# ## For Fibre
#           RMSE        ACC
# [1,] 1.0940558 0.45748208
# [2,] 0.0366813 0.03226115
# > Tab
# PT     RMSE       ACC
# 1  1 1.106064 0.4702890
# 2  2 1.117216 0.4172794
# 3  3 1.045071 0.4683942
# 4  4 1.067748 0.4983423
# 5  5 1.134180 0.4331054
