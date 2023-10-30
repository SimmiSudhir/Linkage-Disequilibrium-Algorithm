##GBLUP using BGLR [5-fold Cross validation]

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
# geno <- as.data.frame(fread("M_dip_1317.txt",header=TRUE))

geno <- as.data.frame(fread("M_con_1317.txt",header=TRUE))

dim(geno) ## 1317 by 25753
geno[1:5,1:5]

clones_1317 <- as.data.frame(fread("1317_clones.txt",header=TRUE))
head(clones_1317)
colnames(clones_1317)[1] <- "Clone"

rownames(geno) <- clones_1317$Clone
geno[1:5,1:5]

## first create a Genomic relation matrix######
## scale the marker matrix
library(BGLR)
X <- scale(geno, center=TRUE, scale=TRUE)
dim(X) ## 1317 25753
p <- ncol(X)
## COMPUTING g
G <- tcrossprod(X)/p
G[1:5,1:5]
dim(G) #  1317 1317


ETA = list(list(K=G, model='RKHS'))

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

## For TCH

# ##RMSE        ACC
# [1,] 8.7975743 0.41101918
# [2,] 0.3162428 0.06354217
# > Tab
# PT     RMSE       ACC
# 1  1 8.825541 0.4357980
# 2  2 8.982996 0.3520779
# 3  3 8.362684 0.5089847
# 4  4 9.183292 0.3644600
# 5  5 8.633359 0.3937753

# ## For CCS
# RMSE       ACC
# [1,] 0.57206382 0.4324571
# [2,] 0.04053885 0.0334592
# > Tab
# PT      RMSE       ACC
# 1  1 0.5895219 0.4089697
# 2  2 0.5063179 0.4510047
# 3  3 0.6139843 0.3860633
# 4  4 0.5660639 0.4656222
# 5  5 0.5844312 0.4506255


# # ## For Fibre (I checked with cont 26 data)
# 
# RMSE        ACC
# [1,] 1.08237338 0.47181318
# [2,] 0.05506656 0.03746291
# > Tab
# PT     RMSE       ACC
# 1  1 1.101274 0.4777194
# 2  2 1.088429 0.5026282
# 3  3 1.042822 0.4852839
# 4  4 1.018732 0.4067901
# 5  5 1.160610 0.4866443

