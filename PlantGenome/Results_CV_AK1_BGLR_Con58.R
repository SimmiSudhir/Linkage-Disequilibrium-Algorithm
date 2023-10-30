###############################################################################
########### RKHS (AK1) ##################################################
############# TCH################################
######## Read Marker Matrix #######################
rm(list=ls())
library(data.table)

geno <- as.data.frame(fread("M_con58_1317.txt",header=TRUE))
dim(geno) ## 1317 by 25753
geno[1:5,1:5]

clones_1317 <- as.data.frame(fread("1317_clones.txt",header=TRUE))
head(clones_1317)
colnames(clones_1317)[1] <- "Clone"

rownames(geno) <- clones_1317$Clone
geno[1:5,1:5]

#######################################################

X <- scale(geno,center=TRUE,scale=TRUE)
p<-ncol(X)

dim(X) ## 1317 25753

x1=X

x2=X

x1tx2 <- x1%*%t(x2)

norm1<-sqrt(apply(x1,1,function(x) crossprod(x)))

norm2<-sqrt(apply(x2,1,function(x) crossprod(x)))

costheta = diag(1/norm1)%*%x1tx2%*%diag(1/norm2)

costheta[which(abs(costheta)>1,arr.ind = TRUE)] = 1

theta<-acos(costheta)

normx1x2<-norm1%*%t(norm2)

J = (sin(theta)+(pi-theta)*cos(theta))

AK1 = 1/pi*normx1x2*J

AK1 <- AK1/median(AK1)

colnames(AK1)<-rownames(x2)

rownames(AK1)<-rownames(x1)
n1<-nrow(AK1)

n2<-ncol(AK1)

nL = 4
AKl1 = AK1

for (l in 1:nL){
  
  AKAK<-tcrossprod(diag(AKl1),diag(AKl1))
  costheta<-AKl1*(AKAK^(-1/2))
  
  costheta[which(costheta>1,arr.ind = TRUE)] = 1
  
  theta<-acos(costheta)
  
  AKl<-(1/pi)*(AKAK^(1/2))*(sin(theta)+(pi-theta)*cos(theta))
  
  AKl1 = AKl
  
}
AKl<-AKl/median(AKl)

rownames(AKl)<-rownames(AK1)

colnames(AKl)<-colnames(AK1)
AKl[1:5,1:5]
dim(AKl) ## 1318 1318

K =AKl
# K =AK1
ETA <- list(K1=list(K=K,model="RKHS"))

###################################################
## Phenotypes ####################

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

## 5 fold cross validation 
library(BGLR)
n <- length(Pheno1$Clone)
K <- 5
set.seed(123)
group <- sample(c(1:K),n, replace=TRUE)
table(group)
#h2_cv <- matrix(data=NA, nrow=K, ncol=1)
Tab <- data.frame(PT=1:K, RMSE=NA, ACC=NA, DIC=NA)

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
  Tab$DIC[g] <- modelBGLR_cv$fit$DIC
}
apply(Tab[-1],2, function(x)c(mean(x),sd(x)))

# ## For TCH
# RMSE        ACC       DIC
# [1,] 8.6440098 0.44470800 7097.7382
# [2,] 0.3450563 0.06544763  157.4809
# > Tab
# PT     RMSE       ACC      DIC
# 1  1 8.609409 0.4806115 6878.295
# 2  2 8.799850 0.3964063 7035.449
# 3  3 8.185555 0.5397896 7122.392
# 4  4 9.115713 0.3804243 7308.475
# 5  5 8.509522 0.4263083 7144.080


# ## For CCS
# 
# RMSE        ACC       DIC
# [1,] 0.56626141 0.45565386 1145.6843
# [2,] 0.04293234 0.04856289  102.6675
# > Tab
# PT      RMSE       ACC      DIC
# 1  1 0.5886225 0.4115430 1031.201
# 2  2 0.4980211 0.4806287 1279.644
# 3  3 0.6068198 0.4079546 1113.275
# 4  4 0.5522358 0.5233877 1222.780
# 5  5 0.5856078 0.4547553 1081.522


# ## For Fibre
# 
#         RMSE        ACC       DIC
# [1,] 1.09270634 0.46269904 2746.2930
# [2,] 0.06522735 0.02568964  131.6041
# > Tab
# PT     RMSE       ACC      DIC
# 1  1 1.126670 0.4518758 2824.233
# 2  2 1.111750 0.4716512 2591.199
# 3  3 1.049776 0.4824694 2830.336
# 4  4 1.005180 0.4229964 2869.463
# 5  5 1.170155 0.4845023 2616.234

######################################################
## for hidden layer 1
rm(list=ls())
library(data.table)

geno <- as.data.frame(fread("M_con58_1317.txt",header=TRUE))
dim(geno) ## 1317 by 25753
geno[1:5,1:5]

clones_1317 <- as.data.frame(fread("1317_clones.txt",header=TRUE))
head(clones_1317)
colnames(clones_1317)[1] <- "Clone"

rownames(geno) <- clones_1317$Clone
geno[1:5,1:5]

#######################################################

X <- scale(geno,center=TRUE,scale=TRUE)
p<-ncol(X)

dim(X) ## 1317 25753

x1=X

x2=X

x1tx2 <- x1%*%t(x2)

norm1<-sqrt(apply(x1,1,function(x) crossprod(x)))

norm2<-sqrt(apply(x2,1,function(x) crossprod(x)))

costheta = diag(1/norm1)%*%x1tx2%*%diag(1/norm2)

costheta[which(abs(costheta)>1,arr.ind = TRUE)] = 1

theta<-acos(costheta)

normx1x2<-norm1%*%t(norm2)

J = (sin(theta)+(pi-theta)*cos(theta))

AK1 = 1/pi*normx1x2*J

AK1 <- AK1/median(AK1)

colnames(AK1)<-rownames(x2)

rownames(AK1)<-rownames(x1)

K =AK1
ETA <- list(K1=list(K=K,model="RKHS"))

###################################################
## Phenotypes ####################

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

## 5 fold cross validation 
library(BGLR)
n <- length(Pheno1$Clone)
K <- 5
set.seed(123)
group <- sample(c(1:K),n, replace=TRUE)
table(group)
#h2_cv <- matrix(data=NA, nrow=K, ncol=1)
Tab <- data.frame(PT=1:K, RMSE=NA, ACC=NA, DIC=NA)

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
  Tab$DIC[g] <- modelBGLR_cv$fit$DIC
}
apply(Tab[-1],2, function(x)c(mean(x),sd(x)))

###########################
# ## For TCH
#       RMSE        ACC       DIC
# [1,] 8.6766203 0.43749810 7409.0296
# [2,] 0.3556138 0.06515341  111.6988
# > Tab
# PT     RMSE       ACC      DIC
# 1  1 8.671859 0.4664819 7264.193
# 2  2 8.864560 0.3842763 7369.627
# 3  3 8.197362 0.5351983 7385.641
# 4  4 9.138492 0.3777066 7564.199
# 5  5 8.510829 0.4238273 7461.488

# ## For CCS
#         RMSE        ACC       DIC
# [1,] 0.56273228 0.46124544 1520.3397
# [2,] 0.04316524 0.04109794   72.0822
# > Tab
# PT      RMSE       ACC      DIC
# 1  1 0.5853105 0.4233191 1453.086
# 2  2 0.4937261 0.4903047 1621.164
# 3  3 0.6033976 0.4197178 1451.500
# 4  4 0.5496285 0.5136362 1557.131
# 5  5 0.5815988 0.4592494 1518.817

# ## For Fibre
# 
#         RMSE        ACC        DIC
# [1,] 1.07441301 0.48468714 2951.38698
# [2,] 0.05651385 0.03408633   45.12731
# > Tab
# PT     RMSE       ACC      DIC
# 1  1 1.102066 0.4807851 2954.309
# 2  2 1.087495 0.5048452 2900.910
# 3  3 1.032578 0.5002879 2977.452
# 4  4 1.003939 0.4270747 3010.501
# 5  5 1.145986 0.5104428 2913.763









