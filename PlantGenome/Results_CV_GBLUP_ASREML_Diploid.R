## fit model
## cross validation[5-FOLDS] GBLUP using asReml  for all traits
## Results are accuracy and RMSE
## diploid model

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

ahatinv<-read.table("GRM_1317_dip.txt",h=T)

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

# ## for TCH
# #         RMSE        ACC
# [1,] 8.8060118 0.40926093
# [2,] 0.3178341 0.06449641
# 
# # > Tab
# PT     RMSE       ACC
# 1  1 8.827269 0.4358816
# 2  2 8.981630 0.3514485
# 3  3 8.382943 0.5075369
# 4  4 9.209118 0.3572810
# 5  5 8.629098 0.3941567
# # 
# # > h2_cv
#         [,1]
# [1,] 0.3757307
# [2,] 0.4256499
# [3,] 0.3336829
# [4,] 0.3841246
# [5,] 0.3998666

# ########## For CCS
#           RMSE        ACC
# # [1,] 0.57219625 0.43212674
# [2,] 0.03989245 0.03184415
# 
# # > Tab
# PT      RMSE       ACC
# 1  1 0.5897921 0.4082111
# 2  2 0.5074076 0.4473001
# 3  3 0.6130571 0.3893773
# 4  4 0.5661757 0.4649718
# 5  5 0.5845487 0.4507734


# # #> h2_cv
#         [,1]
# [1,] 0.5334676
# [2,] 0.5239491
# [3,] 0.5340946
# [4,] 0.5544042
# [5,] 0.5255738

#######################################

## for fibre 
# # #       RMSE        ACC
# [1,] 1.08081999 0.47329289
# [2,] 0.04909754 0.04122598

# # > Tab
# # PT      RMSE       ACC
# 1  1 1.103197 0.4698165
# 2  2 1.079466 0.5146737
# 3  3 1.052336 0.4720653
# 4  4 1.020128 0.4080516
# 5  5 1.148973 0.5018574
# > h2_cv
# # [,1]
# [1,] 0.5566666
# [2,] 0.5823647
# [3,] 0.5574250
# [4,] 0.6261766
# [5,] 0.5687219

##