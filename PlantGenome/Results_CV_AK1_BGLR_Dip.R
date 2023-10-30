
################################################################################
########### RKHS (AK1) ##################################################
############# TCH################################
######## Read Marker Matrix #######################
rm(list=ls())
library(data.table)

geno <- as.data.frame(fread("M_dip_1317.txt",header=TRUE))
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

nL = 3
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


Pheno1 <- cbind(clones, Pheno$TCH)
#Pheno1 <- cbind(clones, Pheno$CCS)
#Pheno1 <- cbind(clones, Pheno$Fibre)

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
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# ##  # for layer 2
#       RMSE        ACC       DIC
# [1,] 8.659368 0.43915474 7181.5582
# [2,] 0.364851 0.07070801  158.2625
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# for layer 3
# RMSE        ACC       DIC
# [1,] 8.670876 0.43772558 7124.8197
# [2,] 0.347729 0.06517131  155.2364

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# ### for layer 6 (TCH)
# 
#         RMSE        ACC       DIC
# [1,] 8.7000706 0.43319827 7010.5997
# [2,] 0.3633444 0.07125243  187.0561

## for layer 7
# #      RMSE       ACC       DIC
# [1,] 8.7147958 0.4310163 6993.9368
# [2,] 0.3487614 0.0673457  167.7835

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# # ## For TCH (for layer 4)
# RMSE       ACC       DIC
# [1,] 8.6767929 0.4373534 7085.0061
# [2,] 0.3496205 0.0678946  173.8627
# > 
#   > Tab
# PT     RMSE       ACC      DIC
# 1  1 8.673933 0.4680784 6842.923
# 2  2 8.833803 0.3849445 6988.880
# 3  3 8.202558 0.5390470 7131.236
# 4  4 9.143554 0.3727963 7295.437
# 5  5 8.530117 0.4219006 7166.554


## For CCS
# ##      RMSE        ACC       DIC
# [1,] 0.5753298 0.42335056 1099.4732
# [2,] 0.0389548 0.04187206   70.6878
# > Tab
# PT      RMSE       ACC      DIC
# 1  1 0.5935261 0.3939056 1044.417
# 2  2 0.5137763 0.4214743 1152.447
# 3  3 0.6139489 0.3823490 1120.897
# 4  4 0.5627440 0.4898939 1171.972
# 5  5 0.5926540 0.4291300 1007.633

# ## For Fibre
# 
#       RMSE        ACC       DIC
# [1,] 1.10307402 0.44305939 2620.6040
# [2,] 0.06308695 0.02725204  206.4729
# > Tab
# PT     RMSE       ACC      DIC
# 1  1 1.128526 0.4357113 2595.352
# 2  2 1.121980 0.4537164 2281.252
# 3  3 1.067956 0.4503394 2770.750
# 4  4 1.015875 0.4010630 2795.607
# 5  5 1.181033 0.4744669 2660.058


###################################################
## for hidden layer 1
rm(list=ls())
library(data.table)

geno <- as.data.frame(fread("M_dip_1317.txt",header=TRUE))
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

# ## For TCH
#       RMSE        ACC        DIC
# [1,] 8.6960479 0.43281832 7401.54849
# [2,] 0.3479174 0.06862898   99.75581
# > Tab
# PT     RMSE       ACC      DIC
# 1  1 8.701689 0.4619109 7280.849
# 2  2 8.899304 0.3719785 7341.812
# 3  3 8.207795 0.5376861 7411.178
# 4  4 9.123329 0.3786322 7546.062
# 5  5 8.548122 0.4138839 7427.840

# ############# For CCS
# 
#         RMSE        ACC        DIC
# [1,] 0.57324467 0.42871736 1448.56401
# [2,] 0.03896263 0.03287968   56.88417
# > 
#   > Tab
# PT      RMSE       ACC      DIC
# 1  1 0.5922234 0.4010145 1384.434
# 2  2 0.5107905 0.4374654 1520.509
# 3  3 0.6121001 0.3919969 1436.733
# 4  4 0.5634879 0.4736146 1492.443
# 5  5 0.5876214 0.4394954 1408.701

# ########## For Fibre
# 
#         RMSE        ACC        DIC
# # [1,] 1.0836429 0.46964530 2839.01814
# # [2,] 0.0498099 0.04236161   93.54023
# 
# Tab
# PT     RMSE       ACC      DIC
# 1  1 1.107579 0.4649827 2829.923
# 2  2 1.091856 0.4990755 2692.143
# 3  3 1.051847 0.4722631 2925.763
# 4  4 1.019152 0.4016004 2915.004
# 5  5 1.147781 0.5103047 2832.257