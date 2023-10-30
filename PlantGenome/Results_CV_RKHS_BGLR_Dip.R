#################################
############ Multi Kernel RKHS regression

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
geno <- as.data.frame(fread("M_dip_1317.txt",header=TRUE))
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

 fwrite(D,file="DM_1317_dip.txt",row.names=F,col.names=T, sep=" ")
## Read the file from the directory ############################
 
 dis_mat <- as.data.frame(fread("DM_1317_dip.txt",header=TRUE))
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
 
# ### For Fibre
#  
#  RMSE        ACC
#  [1,] 1.0954716 0.45697169
#  [2,] 0.0323625 0.02826507
#  > Tab
#  PT     RMSE       ACC
#  1  1 1.098540 0.4828355
#  2  2 1.119499 0.4133978
#  3  3 1.046524 0.4655030
#  4  4 1.084398 0.4775787
#  5  5 1.128397 0.4455436
 
 # ### For TCH
 # RMSE       ACC
 # [1,] 8.7115741 0.4297353
 # [2,] 0.5546087 0.0506280
 # > Tab
 # PT     RMSE       ACC
 # 1  1 8.217657 0.4191384
 # 2  2 8.416487 0.4622260
 # 3  3 9.281347 0.3877854
 # 4  4 8.296470 0.4994454
 # 5  5 9.345909 0.3800810
 
 # ## For CCS
 # RMSE        ACC
 # [1,] 0.56929523 0.44874495
 # [2,] 0.02731051 0.05643073
 # > Tab
 # PT      RMSE       ACC
 # 1  1 0.5683046 0.3883028
 # 2  2 0.5258983 0.4975103
 # 3  3 0.5696024 0.4700195
 # 4  4 0.5833496 0.4996584
 # 5  5 0.5993212 0.3882337