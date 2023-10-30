############ Multi Kernel RKHS regression
## using continuous 58K marker
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
geno <- as.data.frame(fread("M_con58_1317.txt",header=TRUE))
dim(geno) ## 1317 by 58363
geno[1:5,1:5]

clones_1317 <- as.data.frame(fread("1317_clones.txt",header=TRUE))
head(clones_1317)
colnames(clones_1317)[1] <- "Clone"

rownames(geno) <- clones_1317$Clone
geno[1:5,1:5]

## scale the marker matrix
library(BGLR)
X <- scale(geno, center=TRUE, scale=TRUE)
dim(X) ## 1317 58363
p <- ncol(X)

## calculate euclidean matrix #######

D <- (as.matrix(dist(X,method='euclidean'))^2)/p
D[1:5,1:5]
dim(D)
######## write the distance matrix in directory

fwrite(D,file="DM_1317_con58.txt",row.names=F,col.names=T, sep=" ")
## Read the file from the directory ############################

dis_mat <- as.data.frame(fread("DM_1317_con58.txt",header=TRUE))
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

# ## For TCH
# RMSE        ACC
# [1,] 8.6999844 0.43183861
# [2,] 0.6214385 0.05801931
# > Tab
# PT     RMSE       ACC
# 1  1 8.096340 0.4477954
# 2  2 8.398038 0.4676930
# 3  3 9.281912 0.3877368
# 4  4 8.270983 0.4985180
# 5  5 9.452649 0.3574498

# ## For CCS
# RMSE        ACC
# [1,] 0.56243814 0.47156054
# [2,] 0.02295617 0.05182418
# > Tab
# PT      RMSE       ACC
# 1  1 0.5591235 0.4187322
# 2  2 0.5277127 0.4921935
# 3  3 0.5631454 0.4909945
# 4  4 0.5712627 0.5373858
# 5  5 0.5909465 0.4184966

# ## For Fibre
# 
#         RMSE        ACC
# [1,] 1.08346504 0.47544331
# [2,] 0.03151431 0.02571104
# > Tab
# PT     RMSE       ACC
# 1  1 1.096148 0.4849752
# 2  2 1.106111 0.4363898
# 3  3 1.039597 0.4764647
# 4  4 1.061868 0.5073148
# 5  5 1.113600 0.4720720







