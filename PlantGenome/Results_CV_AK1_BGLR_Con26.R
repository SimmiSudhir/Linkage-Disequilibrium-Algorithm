###############################################################################
########### RKHS (AK1) ##################################################
############# TCH################################
######## Read Marker Matrix #######################
rm(list=ls())
library(data.table)

geno <- as.data.frame(fread("M_con_1317.txt",header=TRUE))
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

###
# ## For TCH
# 
# RMSE        ACC       DIC
# [1,] 8.6736511 0.43875142 7139.7484
# [2,] 0.3515836 0.06750187  144.0118
# > Tab
# PT     RMSE       ACC      DIC
# 1  1 8.669721 0.4684444 6943.001
# 2  2 8.803173 0.3956525 7079.030
# 3  3 8.196224 0.5410019 7162.283
# 4  4 9.156017 0.3710376 7336.136
# 5  5 8.543120 0.4176207 7178.292

# ## For CCS
# 
# RMSE        ACC       DIC
# [1,] 0.56808796 0.44856986 1195.3696
# [2,] 0.04469276 0.05347616  125.1355
# > Tab
# PT      RMSE       ACC      DIC
# 1  1 0.5923043 0.3981051 1058.802
# 2  2 0.4968687 0.4855323 1361.444
# 3  3 0.6077163 0.4044539 1182.550
# 4  4 0.5526237 0.5216494 1276.749
# 5  5 0.5909267 0.4331085 1097.303

# ## For Fibre
# 
# RMSE        ACC       DIC
# [1,] 1.09913799 0.45170365 2790.0980
# [2,] 0.06667968 0.02716234  115.7789
# > Tab
# PT     RMSE       ACC      DIC
# 1  1 1.128218 0.4492859 2874.611
# 2  2 1.115649 0.4663402 2676.663
# 3  3 1.054255 0.4750048 2860.455
# 4  4 1.013193 0.4060392 2887.280
# 5  5 1.184375 0.4618481 2651.481


######################################################
## for hidden layer 1
rm(list=ls())
library(data.table)

geno <- as.data.frame(fread("M_con_1317.txt",header=TRUE))
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

# ### For TCH
# 
#         RMSE       ACC       DIC
# [1,] 8.7139341 0.4295998 7442.6212
# # [2,] 0.3488486 0.0666583  102.1874
# 
# 
# PT     RMSE       ACC      DIC
# 1  1 8.706853 0.4598504 7307.790
# 2  2 8.892576 0.3744993 7402.683
# 3  3 8.249465 0.5294757 7432.267
# 4  4 9.173589 0.3693656 7583.927
# 5  5 8.547188 0.4148078 7486.438


# ############ For CCS
# 
#         RMSE        ACC        DIC
# [1,] 0.56537693 0.45336225 1558.75818
# [2,] 0.04537553 0.04807561   79.25999
# > Tab
# PT      RMSE       ACC      DIC
# 1  1 0.5909086 0.4050878 1484.328
# 2  2 0.4929510 0.4931806 1668.034
# 3  3 0.6061052 0.4118240 1483.982
# 4  4 0.5502862 0.5126668 1603.812
# 5  5 0.5866337 0.4440520 1553.636

# ## For Fibre
#         RMSE        ACC       DIC
# [1,] 1.08275489 0.47196371 2993.8456
# [2,] 0.05972001 0.03325891   37.3625
# > Tab
# PT     RMSE       ACC      DIC
# 1  1 1.106564 0.4748024 3002.442
# 2  2 1.094490 0.4950609 2964.649
# 3  3 1.038531 0.4919915 3005.745
# 4  4 1.011042 0.4141563 3045.514
# 5  5 1.163147 0.4838075 2950.878







