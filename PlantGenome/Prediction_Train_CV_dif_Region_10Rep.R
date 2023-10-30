rm(list=ls())

# evaluate the model's performance separately for each region

#######################################################
library(data.table)
data <- as.data.frame(fread("pheno_1318.txt", header=TRUE))
str(data)
data$Series <- as.factor(data$Series)
data$Region <- as.factor(data$Region)
data$Trial <- as.factor(data$Trial)
data$Crop <- as.factor(data$Crop)
data$Clone <- as.factor(data$Clone)
str(data)
### exclude 2013 year data set 
Data <- subset(data, Series!=2013)
Data <- droplevels(Data)
str(Data) 

############
library(caret)

# Get the unique regions in the dataset
unique_regions <- unique(Data$Region)
#########################################################
################ Region A ##############################

# Subset the data for the current region
region_data <- subset(Data, Region == "C")
region_data <- droplevels(region_data)
str(region_data)

#######################################################
library(data.table)
Z <- as.data.frame(fread("M_dip_1317.txt",header=TRUE))
dim(Z) ## 1317, 25753
Z[1:5,1:5]
clones <- read.table ("1317_clones.txt", header = T)
head(clones)
colnames(clones)[1] <- "clone"

rownames(Z) <- clones$clone

Z_201_cl <- intersect(rownames(Z),region_data$Clone)
length(Z_201_cl) ## 201

Z_201 <- Z[rownames(Z)%in% Z_201_cl ,]
dim(Z_201) ## 201 25753
Z_201[1:5,1:5]

any(is.na(Z_201)) ## FALSE
equal_cols <- which(apply(Z_201, 2, function(x) sd(x) == 0))
Z_201_fil <- Z_201[,-equal_cols]

dim(Z_201_fil) ## 201 25748
X_201 <- scale(Z_201_fil, center=TRUE, scale=TRUE)
# any(is.na(X_379)) ## TRUE
any(is.na(X_201)) ## FALSE

#X_201 <- scale(Z_201, center=TRUE, scale=TRUE)
dim(X_201) ## 201 25748

p <- ncol(X_201)
## COMPUTING g
G_201 <- tcrossprod(X_201)/p
G_201[1:5,1:5]
dim(G_201) ## 201 201

#######################################################
det(G_201)# -1.161617e-37
for (i in 1:100) {  # i added 8x 0.05 in total
  
  diag(G_201)=diag(G_201)+0.05
  d = det(G_201)
  if (d!=0) break 
  
}
print (i) ####### 1
print (d)## 2.824176e-19

########################################################

########################################################################################
Ginv <- solve(G_201)
Ginv[1:5,1:5]
dim(Ginv) # 201 201
###################################################################
source('C:/Users/uqsyada1/OneDrive - The University of Queensland/PhD/Trainings and workshop/201910_AsREML/Asreml-R workshop/Asreml(4-5)/Material_UQ45/full2sparse.R')
Ginv.sparse<-full2sparse(Ginv)
colnames(Ginv.sparse)<-c('Row','Column','Ainverse')
head(Ginv.sparse)

##################################################

class(Z_201_cl) ## character
cl_201 <- as.data.frame(Z_201_cl)
colnames(cl_201)[1]<-"clone"



######################################################
ahatinv <- Ginv.sparse
ahatinv<-as.matrix.data.frame(ahatinv)
attr(ahatinv,"rowNames")<-as.character(cl_201$clone)
attr(ahatinv,"colNames")<-as.character(cl_201$clone)
attr(ahatinv,"INVERSE")<-TRUE
head(ahatinv)



######################################################
library(asreml)
# Set a seed for reproducibility
set.seed(42)
# Initialize a vector to store accuracy results
accuracy_results <- numeric(10)
for (iter in 1:10) {
  # Get the unique clones
  unique_clones <- unique(region_data$Clone)
  # Generate random indices for the test set
  test_indices <- sample(length(unique_clones), floor(0.2 * length(unique_clones)))
  # Create the test set by subsetting the data based on the selected clones and setting the response variable to NA
  test_set <- region_data[region_data$Clone %in% unique_clones[test_indices], ]
  test_set <- droplevels(test_set)
  str(test_set) ## 113 clones
  
  train_data <- region_data
  exclude <- intersect(train_data$Clone, test_set$Clone) ## 113
  train_data$TCHBlup [train_data$Clone %in% exclude] = NA
  train_data$CCSBlup [train_data$Clone %in% exclude] = NA
  train_data$FibreBlup [train_data$Clone %in% exclude] = NA
##################### TCH ######################################################

  modelGBLUP_A = asreml(fixed=FibreBlup ~ Series + Region + Trial + Crop,
                        random=~vm(Clone,ahatinv),
                        #random=~Clone,
                        workspace=128e06,na.action = na.method(y="include"),
                        data= train_data)
  summary(modelGBLUP_A)$varcomp
  source('C:/Users/uqsyada1/OneDrive - The University of Queensland/PhD/Trainings and workshop/201910_AsREML/Asreml-R workshop/ASReml(4-5)/Material_UQ45/pin.r')
  pin(modelGBLUP_A,h2_1~V1/(V1+V2))
  ##########################################################################
  
  ##############################################
  BLUPS_Fibre <- summary(modelGBLUP_A,coef=TRUE)$coef.random
  BLUPS_Fibre = as.data.frame(BLUPS_Fibre)
  head(BLUPS_Fibre)
  dim(BLUPS_Fibre) # 201    3
  a.effects <- BLUPS_Fibre[1:201,1]
  #######################################################
  
  # getting prediction accuracies, i.e. correlating actual with predicted performance
  
  all.names = c()
  for (i in 1:nrow (BLUPS_Fibre) ) {
    
    name = strsplit (rownames(BLUPS_Fibre[i,]), split = "_") [[1]][2]
    
    all.names = c(all.names, name)
    
  }
  
  BLUPS_Fibre$clone = all.names
  head(BLUPS_Fibre)
  
  test_set$GEBV_Fibre = 0
  
  for (i in 1:nrow(test_set)){
    
    clone = noquote(test_set[i,'Clone'])
    
    sum <-BLUPS_Fibre$solution[BLUPS_Fibre$clone == clone][1] 
    
    test_set$GEBV_Fibre[test_set$Clone == clone] = sum
    
  }
  ###################################################################################
  head(test_set,5)
  accuracy <- cor (test_set$FibreBlup, test_set$GEBV_Fibre,use = "pairwise.complete.obs")## 0.3235673, 0.1469458, 0.2054209
 
accuracy_results[iter] <- accuracy
}
################# Region A (TCH)
# ## accuracy_results
# [1] 0.3111043 0.4218889 0.3797670 0.2497096 0.3988971
# [6] 0.2602388 0.3147351 0.3289753 0.2550326 0.2516247
mean(accuracy_results) # 0.3171973
sd(accuracy_results)## 0.06471746
# Calculate standard error of accuracy results
se_accuracy <- sd(accuracy_results) / sqrt(length(accuracy_results)) ## 0.02046546
#########################################################

##############################################################
############## ################# Region A (CCS)
# #accuracy_results
# [1] 0.3534350 0.1765214 0.2191342 0.1729970 0.1188205
# [6] 0.1710285 0.2632158 0.1681101 0.2734341 0.1468880
mean(accuracy_results) # 0.2063585
sd(accuracy_results)## 0.07105602
# Calculate standard error of accuracy results
se_accuracy <- sd(accuracy_results) / sqrt(length(accuracy_results)) ## 0.02246989
#########################################################

# ###################### Region A (Fibre)
# #accuracy_results
# [1] 0.3258270 0.3854530 0.3072652 0.4484732 0.3839438
# [6] 0.4630164 0.4594346 0.3758794 0.3485980 0.3759429
mean(accuracy_results) ## 0.39
sd(accuracy_results) ## 0.05
se_accuracy <- sd(accuracy_results) / sqrt(length(accuracy_results))## 0.017

##############################################################
###############################################################
########### Region S (TCH)
# ####accuracy_results
# [1] 0.3240686 0.3142145 0.3147633 0.2669639 0.4719881
# [6] 0.2793157 0.3863270 0.2950498 0.2593800 0.3253976
mean(accuracy_results) ## 0.32
sd(accuracy_results) ## 0.06
se_accuracy <- sd(accuracy_results) / sqrt(length(accuracy_results))## 0.020

# ########### Region S (CCS)
# accuracy_results
# [1] 0.3440145 0.4411625 0.4293758 0.3674983 0.4467392
# [6] 0.2643457 0.4389816 0.4438984 0.4197657 0.3206335
mean(accuracy_results) ## 0.39
sd(accuracy_results) ## 0.06
se_accuracy <- sd(accuracy_results) / sqrt(length(accuracy_results))## 0.020

# # ########### Region S (Fibre)
# accuracy_results
# [1] 0.4280605 0.3320165 0.3205004 0.2075396 0.1927608
# [6] 0.1152725 0.1962196 0.3517864 0.2406572 0.2535757
mean(accuracy_results) ## 0.26
sd(accuracy_results) ## 0.09
se_accuracy <- sd(accuracy_results) / sqrt(length(accuracy_results))## 0.029

####################################################
# ################ Region N (TCH)
# accuracy_results
# [1] 0.14212103 0.18646280 0.18216613 0.24168409
# [5] 0.23868671 0.06276132 0.33720793 0.17094191
# [9] 0.31042855 0.14491779
mean(accuracy_results) ## 0.20
sd(accuracy_results) ## 0.08
se_accuracy <- sd(accuracy_results) / sqrt(length(accuracy_results))## 0.026
range(accuracy_results) ## (0.063, 0.337)

################ Region N (CCS)
# accuracy_results
# [1] 0.2233292 0.4386351 0.3596124 0.3738894 0.4362202
#[6] 0.2827107 0.3746723 0.3101184 0.4643429 0.2268781
mean(accuracy_results) ## 0.35
sd(accuracy_results) ## 0.09
se_accuracy <- sd(accuracy_results) / sqrt(length(accuracy_results))## 0.027
range(accuracy_results) ## (0.2233292 0.4643429)

# ################ Region N (Fibre)
# accuracy_results
# 
# [1] 0.3989093 0.5459054 0.4190848 0.2145253 0.4421127
# [6] 0.3569853 0.4464605 0.4408533 0.3517553 0.4808933
mean(accuracy_results) ## 0.41
sd(accuracy_results) ## 0.09
se_accuracy <- sd(accuracy_results) / sqrt(length(accuracy_results))## 0.028
range(accuracy_results) ## (0.2145253 0.5459054)

##########################################
##########################################
# ### Region C (TCH)
# accuracy_results
# [1]  0.13815765 -0.01463427  0.23332882  0.15031822
# [5] -0.05565870  0.04465349 -0.01440036  0.24888747
# [9]  0.19717485  0.22865784
mean(accuracy_results) ## 0.1156485
sd(accuracy_results) ## 0.1159503
se_accuracy <- sd(accuracy_results) / sqrt(length(accuracy_results))## 0.037
range(accuracy_results) # (-0.0556587  0.2488875)

############ Region C (CCS)
accuracy_results
#[1] 0.2217085 0.3273817 0.5579216 0.3677489 0.3500606
#[6] 0.3119095 0.1876895 0.3281250 0.4540205 0.4679703
mean(accuracy_results) ## 0.36
sd(accuracy_results) ## 0.11
se_accuracy <- sd(accuracy_results) / sqrt(length(accuracy_results))## 0.035
range(accuracy_results) # (0.1876895 0.5579216)

# ############ Region C (Fibre)
# accuracy_results
# [1] 0.14694584 0.20542088 0.28029039 0.07988520
# [5] 0.38073256 0.25001131 0.24379859 0.32525916
# [9] 0.07400035 0.25578689
mean(accuracy_results) ## 0.22
sd(accuracy_results) ## 0.0998
se_accuracy <- sd(accuracy_results) / sqrt(length(accuracy_results))## 0.032
range(accuracy_results) # (0.07400035 0.38073256)

