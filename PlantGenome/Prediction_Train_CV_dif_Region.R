
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

  # Subset the data for the current region
  region_data <- subset(Data, Region == "A")
  region_data <- droplevels(region_data)
  str(region_data)
  # Set a seed for reproducibility
  set.seed(42)
  
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

#######################################################
library(data.table)
Z <- as.data.frame(fread("M_dip_1317.txt",header=TRUE))
dim(Z) ## 1317, 25753
Z[1:5,1:5]
clones <- read.table ("1317_clones.txt", header = T)
head(clones)
colnames(clones)[1] <- "clone"

rownames(Z) <- clones$clone

Z_567_cl <- intersect(rownames(Z),region_data$Clone)
length(Z_567_cl) ## 567

Z_567 <- Z[rownames(Z)%in% Z_567_cl ,]
dim(Z_567) ## 567 25753
Z_567[1:5,1:5]

X_567 <- scale(Z_567, center=TRUE, scale=TRUE)
dim(X_567) ## 540 25753
p <- ncol(X_567)
## COMPUTING g
G_567 <- tcrossprod(X_567)/p
G_567[1:5,1:5]
dim(G_567) ## 567 567

#######################################################
det(G_567)# 2.560583e-153
for (i in 1:100) {  # i added 8x 0.05 in total
  
  diag(G_567)=diag(G_567)+0.05
  d = det(G_567)
  if (d!=0) break 
  
}
print (i) ####### 1
print (d)## 2.499913e-112

########################################################

########################################################################################
Ginv <- solve(G_567)
Ginv[1:5,1:5]
dim(Ginv) # 567 567
###################################################################
source('C:/Users/uqsyada1/OneDrive - The University of Queensland/PhD/Trainings and workshop/201910_AsREML/Asreml-R workshop/Asreml(4-5)/Material_UQ45/full2sparse.R')
Ginv.sparse<-full2sparse(Ginv)
colnames(Ginv.sparse)<-c('Row','Column','Ainverse')
head(Ginv.sparse)

##################################################

class(Z_567_cl) ## character
cl_567 <- as.data.frame(Z_567_cl)
colnames(cl_567)[1]<-"clone"



######################################################
ahatinv <- Ginv.sparse
ahatinv<-as.matrix.data.frame(ahatinv)
attr(ahatinv,"rowNames")<-as.character(cl_567$clone)
attr(ahatinv,"colNames")<-as.character(cl_567$clone)
attr(ahatinv,"INVERSE")<-TRUE
head(ahatinv)

######################################################
library(asreml)
##################### TCH ######################################################


modelGBLUP_A = asreml(fixed=TCHBlup ~ Series + Region + Trial + Crop,
                      random=~vm(Clone,ahatinv),
                      #random=~Clone,
                      workspace=128e06,na.action = na.method(y="include"),
                      data= train_data)
summary(modelGBLUP_A)$varcomp
# #
# component std.error  z.ratio bound %ch
# vm(Clone, ahatinv) 132.13635 9.3907813 14.07086     P   0
# units!R             29.92912 0.7491026 39.95330     P   0

source('C:/Users/uqsyada1/OneDrive - The University of Queensland/PhD/Trainings and workshop/201910_AsREML/Asreml-R workshop/ASReml(4-5)/Material_UQ45/pin.r')
pin(modelGBLUP_A,h2_1~V1/(V1+V2))
# #
# Estimate        SE
# h2_1 0.815327 0.0114436
##########################################################################

##############################################
BLUPS_TCH <- summary(modelGBLUP_A,coef=TRUE)$coef.random
BLUPS_TCH = as.data.frame(BLUPS_TCH)
head(BLUPS_TCH)
dim(BLUPS_TCH) # 567    3
a.effects <- BLUPS_TCH[1:567,1]
#######################################################

# getting prediction accuracies, i.e. correlating actual with predicted performance

all.names = c()
for (i in 1:nrow (BLUPS_TCH) ) {
  
  name = strsplit (rownames(BLUPS_TCH[i,]), split = "_") [[1]][2]
  
  all.names = c(all.names, name)
  
}

BLUPS_TCH$clone = all.names
head(BLUPS_TCH)

test_set$GEBV_TCH = 0

for (i in 1:nrow(test_set)){
  
  clone = noquote(test_set[i,'Clone'])
  
  sum <-BLUPS_TCH$solution[BLUPS_TCH$clone == clone][1] 
  
  test_set$GEBV_TCH[test_set$Clone == clone] = sum
  
}
###################################################################################
head(test_set,5)
cor (test_set$TCHBlup, test_set$GEBV_TCH)## 0.3111043
plot(test_set$TCHBlup, test_set$GEBV_TCH)

############################################################
## CCS

exclude <- intersect(train_data$Clone, test_set$Clone) ## 113
train_data$CCSBlup [train_data$Clone %in% exclude] = NA


modelGBLUP_A = asreml(fixed=CCSBlup ~ Series + Region + Trial + Crop,
                      random=~vm(Clone,ahatinv),
                      #random=~Clone,
                      workspace=128e06,na.action = na.method(y="include"),
                      data= train_data)
summary(modelGBLUP_A)$varcomp

# ##
# component   std.error  z.ratio bound %ch
# vm(Clone, ahatinv) 0.55582391 0.038669084 14.37386     P   0
# units!R            0.08793311 0.002199435 39.97986     P   0

source('C:/Users/uqsyada1/OneDrive - The University of Queensland/PhD/Trainings and workshop/201910_AsREML/Asreml-R workshop/ASReml(4-5)/Material_UQ45/pin.r')
pin(modelGBLUP_A,h2_1~V1/(V1+V2))
# #
# Estimate          SE
# h2_1 0.8634064 0.008770518
##########################################################################

##############################################
BLUPS_CCS <- summary(modelGBLUP_A,coef=TRUE)$coef.random
BLUPS_CCS = as.data.frame(BLUPS_CCS)
head(BLUPS_CCS)
dim(BLUPS_CCS) # 567    3
a.effects <- BLUPS_CCS[1:567,1]
#######################################################

# getting prediction accuracies, i.e. correlating actual with predicted performance

all.names = c()
for (i in 1:nrow (BLUPS_CCS) ) {
  
  name = strsplit (rownames(BLUPS_CCS[i,]), split = "_") [[1]][2]
  
  all.names = c(all.names, name)
  
}

BLUPS_CCS$clone = all.names
head(BLUPS_CCS)

test_set$GEBV_CCS = 0

for (i in 1:nrow(test_set)){
  
  clone = noquote(test_set[i,'Clone'])
  
  sum <-BLUPS_CCS$solution[BLUPS_CCS$clone == clone][1] 
  
  test_set$GEBV_CCS[test_set$Clone == clone] = sum
  
}
###################################################################################
head(test_set,5)
cor (test_set$CCSBlup, test_set$GEBV_CCS)##  0.353435
plot(test_set$CCSBlup, test_set$GEBV_CCS)

############################################################
### Fibre 

exclude <- intersect(train_data$Clone, test_set$Clone) ## 113
train_data$FibreBlup [train_data$Clone %in% exclude] = NA


modelGBLUP_A = asreml(fixed=FibreBlup ~ Series + Region + Trial + Crop,
                      random=~vm(Clone,ahatinv),
                      #random=~Clone,
                      workspace=128e06,na.action = na.method(y="include"),
                      data= train_data)
summary(modelGBLUP_A)$varcomp

# ##
# component   std.error  z.ratio bound %ch
# vm(Clone, ahatinv) 2.1976918 0.148350416 14.81419     P   0
# units!R            0.1018235 0.002546356 39.98791     P   0


source('C:/Users/uqsyada1/OneDrive - The University of Queensland/PhD/Trainings and workshop/201910_AsREML/Asreml-R workshop/ASReml(4-5)/Material_UQ45/pin.r')
pin(modelGBLUP_A,h2_1~V1/(V1+V2))

# #
# Estimate         SE
# h2_1 0.9557196 0.00305347
##########################################################################

##############################################
BLUPS_Fibre <- summary(modelGBLUP_A,coef=TRUE)$coef.random
BLUPS_Fibre = as.data.frame(BLUPS_Fibre)
head(BLUPS_Fibre)
dim(BLUPS_Fibre) # 567    3
a.effects <- BLUPS_Fibre[1:567,1]
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
cor (test_set$FibreBlup, test_set$GEBV_Fibre)## 0.325827
plot(test_set$FibreBlup, test_set$GEBV_Fibre)
#######################################################
##########################################################
### Region S
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Subset the data for the current region
region_data <- subset(Data, Region == "S")
region_data <- droplevels(region_data)
str(region_data) ## 411 LEVELS
# Set a seed for reproducibility
set.seed(42)

# Get the unique clones
unique_clones <- unique(region_data$Clone)
# Generate random indices for the test set
test_indices <- sample(length(unique_clones), floor(0.2 * length(unique_clones)))
# Create the test set by subsetting the data based on the selected clones and setting the response variable to NA
test_set <- region_data[region_data$Clone %in% unique_clones[test_indices], ]
test_set <- droplevels(test_set)
str(test_set) ## 82 clones

train_data <- region_data
exclude <- intersect(train_data$Clone, test_set$Clone) ## 82
train_data$TCHBlup [train_data$Clone %in% exclude] = NA

#######################################################
library(data.table)
Z <- as.data.frame(fread("M_dip_1317.txt",header=TRUE))
dim(Z) ## 1317, 25753
Z[1:5,1:5]
clones <- read.table ("1317_clones.txt", header = T)
head(clones)
colnames(clones)[1] <- "clone"

rownames(Z) <- clones$clone

Z_411_cl <- intersect(rownames(Z),region_data$Clone)
length(Z_411_cl) ## 411

Z_411 <- Z[rownames(Z)%in% Z_411_cl ,]
dim(Z_411) ## 411 25753
Z_411[1:5,1:5]

X_411 <- scale(Z_411, center=TRUE, scale=TRUE)
dim(X_411) ## 411 25753
p <- ncol(X_411)
## COMPUTING g
G_411 <- tcrossprod(X_411)/p
G_411[1:5,1:5]
dim(G_411) ## 411 411

#######################################################
det(G_411)# -5.304497e-93
for (i in 1:100) {  # i added 8x 0.05 in total
  
  diag(G_411)=diag(G_411)+0.05
  d = det(G_411)
  if (d!=0) break 
  
}
print (i) ####### 1
print (d)## 1.523998e-63

########################################################

########################################################################################
Ginv <- solve(G_411)
Ginv[1:5,1:5]
dim(Ginv) # 567 567
###################################################################
source('C:/Users/uqsyada1/OneDrive - The University of Queensland/PhD/Trainings and workshop/201910_AsREML/Asreml-R workshop/Asreml(4-5)/Material_UQ45/full2sparse.R')
Ginv.sparse<-full2sparse(Ginv)
colnames(Ginv.sparse)<-c('Row','Column','Ainverse')
head(Ginv.sparse)

##################################################

class(Z_411_cl) ## character
cl_411 <- as.data.frame(Z_411_cl)
colnames(cl_411)[1]<-"clone"



######################################################
ahatinv <- Ginv.sparse
ahatinv<-as.matrix.data.frame(ahatinv)
attr(ahatinv,"rowNames")<-as.character(cl_411$clone)
attr(ahatinv,"colNames")<-as.character(cl_411$clone)
attr(ahatinv,"INVERSE")<-TRUE
head(ahatinv)

######################################################
library(asreml)
##################### TCH ######################################################


modelGBLUP_A = asreml(fixed=TCHBlup ~ Series + Region + Trial + Crop,
                      random=~vm(Clone,ahatinv),
                      #random=~Clone,
                      workspace=128e06,na.action = na.method(y="include"),
                      data= train_data)
summary(modelGBLUP_A)$varcomp

# # #
# #    component std.error  z.ratio bound %ch
# vm(Clone, ahatinv) 116.48354 9.4973598 12.26483     P   0
# units!R             36.95564 0.7851202 47.07004     P   0

source('C:/Users/uqsyada1/OneDrive - The University of Queensland/PhD/Trainings and workshop/201910_AsREML/Asreml-R workshop/ASReml(4-5)/Material_UQ45/pin.r')
pin(modelGBLUP_A,h2_1~V1/(V1+V2))

# #Estimate         SE
# h2_1 0.7591512 0.01544986
##########################################################################

##############################################
BLUPS_TCH <- summary(modelGBLUP_A,coef=TRUE)$coef.random
BLUPS_TCH = as.data.frame(BLUPS_TCH)
head(BLUPS_TCH)
dim(BLUPS_TCH) # 411    3
a.effects <- BLUPS_TCH[1:411,1]
#######################################################

# getting prediction accuracies, i.e. correlating actual with predicted performance

all.names = c()
for (i in 1:nrow (BLUPS_TCH) ) {
  
  name = strsplit (rownames(BLUPS_TCH[i,]), split = "_") [[1]][2]
  
  all.names = c(all.names, name)
  
}

BLUPS_TCH$clone = all.names
head(BLUPS_TCH)

test_set$GEBV_TCH = 0

for (i in 1:nrow(test_set)){
  
  clone = noquote(test_set[i,'Clone'])
  
  sum <-BLUPS_TCH$solution[BLUPS_TCH$clone == clone][1] 
  
  test_set$GEBV_TCH[test_set$Clone == clone] = sum
  
}
###################################################################################
head(test_set,5)
cor (test_set$TCHBlup, test_set$GEBV_TCH,use = "pairwise.complete.obs")## 0.3249668
plot(test_set$TCHBlup, test_set$GEBV_TCH)

############################################################
## CCS

exclude <- intersect(train_data$Clone, test_set$Clone) ## 113
train_data$CCSBlup [train_data$Clone %in% exclude] = NA


modelGBLUP_A = asreml(fixed=CCSBlup ~ Series + Region + Trial + Crop,
                      random=~vm(Clone,ahatinv),
                      #random=~Clone,
                      workspace=128e06,na.action = na.method(y="include"),
                      data= train_data)
summary(modelGBLUP_A)$varcomp
# # #
# # c component   std.error  z.ratio bound %ch
# vm(Clone, ahatinv) 0.5645910 0.045456522 12.42046     P   0
# units!R            0.1365703 0.002848088 47.95157     P   0  0

source('C:/Users/uqsyada1/OneDrive - The University of Queensland/PhD/Trainings and workshop/201910_AsREML/Asreml-R workshop/ASReml(4-5)/Material_UQ45/pin.r')
pin(modelGBLUP_A,h2_1~V1/(V1+V2))
# #Estimate        SE
# h2_10.8052227 0.01307
##########################################################################

##############################################
BLUPS_CCS <- summary(modelGBLUP_A,coef=TRUE)$coef.random
BLUPS_CCS = as.data.frame(BLUPS_CCS)
head(BLUPS_CCS)
dim(BLUPS_CCS) # 411    3
a.effects <- BLUPS_CCS[1:411,1]
#######################################################

# getting prediction accuracies, i.e. correlating actual with predicted performance

all.names = c()
for (i in 1:nrow (BLUPS_CCS) ) {
  
  name = strsplit (rownames(BLUPS_CCS[i,]), split = "_") [[1]][2]
  
  all.names = c(all.names, name)
  
}

BLUPS_CCS$clone = all.names
head(BLUPS_CCS)

test_set$GEBV_CCS = 0

for (i in 1:nrow(test_set)){
  
  clone = noquote(test_set[i,'Clone'])
  
  sum <-BLUPS_CCS$solution[BLUPS_CCS$clone == clone][1] 
  
  test_set$GEBV_CCS[test_set$Clone == clone] = sum
  
}
###################################################################################
head(test_set,5)
cor (test_set$CCSBlup, test_set$GEBV_CCS)##0.3440145
plot(test_set$CCSBlup, test_set$GEBV_CCS)

############################################################
### Fibre 

exclude <- intersect(train_data$Clone, test_set$Clone) ## 113
train_data$FibreBlup [train_data$Clone %in% exclude] = NA


modelGBLUP_A = asreml(fixed=FibreBlup ~ Series + Region + Trial + Crop,
                      random=~vm(Clone,ahatinv),
                      #random=~Clone,
                      workspace=128e06,na.action = na.method(y="include"),
                      data= train_data)
summary(modelGBLUP_A)$varcomp
# # #
# # #   component   std.error  z.ratio bound %ch
# vm(Clone, ahatinv) 1.29575473 0.102316639 12.66416     P   0
# units!R            0.09916748 0.002153174 46.05641     P   0

source('C:/Users/uqsyada1/OneDrive - The University of Queensland/PhD/Trainings and workshop/201910_AsREML/Asreml-R workshop/ASReml(4-5)/Material_UQ45/pin.r')
pin(modelGBLUP_A,h2_1~V1/(V1+V2))
# #
# Estimate          SE
# h2_1 0.9289082 0.005412696
##########################################################################

##############################################
BLUPS_Fibre <- summary(modelGBLUP_A,coef=TRUE)$coef.random
BLUPS_Fibre = as.data.frame(BLUPS_Fibre)
head(BLUPS_Fibre)
dim(BLUPS_Fibre) # 411   3
a.effects <- BLUPS_Fibre[1:411,1]
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
cor (test_set$FibreBlup, test_set$GEBV_Fibre,use = "pairwise.complete.obs")## 0.4280605
plot(test_set$FibreBlup, test_set$GEBV_Fibre)

##########################################################
### Region N
##################################################
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Subset the data for the current region
region_data <- subset(Data, Region == "N")
region_data <- droplevels(region_data)
str(region_data) ## 379 LEVELS
# Set a seed for reproducibility
set.seed(42)

# Get the unique clones
unique_clones <- unique(region_data$Clone)
# Generate random indices for the test set
test_indices <- sample(length(unique_clones), floor(0.2 * length(unique_clones)))
# Create the test set by subsetting the data based on the selected clones and setting the response variable to NA
test_set <- region_data[region_data$Clone %in% unique_clones[test_indices], ]
test_set <- droplevels(test_set)
str(test_set) ## 75 clones

train_data <- region_data
exclude <- intersect(train_data$Clone, test_set$Clone) ## 75
train_data$TCHBlup [train_data$Clone %in% exclude] = NA

#######################################################
library(data.table)
Z <- as.data.frame(fread("M_dip_1317.txt",header=TRUE))
dim(Z) ## 1317, 25753
Z[1:5,1:5]
clones <- read.table ("1317_clones.txt", header = T)
head(clones)
colnames(clones)[1] <- "clone"

rownames(Z) <- clones$clone

Z_379_cl <- intersect(rownames(Z),region_data$Clone)
length(Z_379_cl) ##379

Z_379 <- Z[rownames(Z)%in% Z_379_cl ,]
dim(Z_379) ## 379 25753
Z_379[1:5,1:5]
any(is.na(Z_379)) ## FALSE
equal_cols <- which(apply(Z_379, 2, function(x) sd(x) == 0))
Z_379_fil <- Z_379[,-equal_cols]

dim(Z_379_fil) ## 379 25748
X_379 <- scale(Z_379_fil, center=TRUE, scale=TRUE)
# any(is.na(X_379)) ## TRUE
any(is.na(X_379)) ## FALSE

# Find the NA values in X_379
#na_indices <- which(is.na(X_379), arr.ind = TRUE)
dim(X_379) ## 379 25748
p <- ncol(X_379)
## COMPUTING g
G_379 <- tcrossprod(X_379)/p
G_379[1:5,1:5]
dim(G_379) ## 379 379

#######################################################
det(G_379)# -1.66497e-80
for (i in 1:100) {  # i added 8x 0.05 in total
  
  diag(G_379)=diag(G_379)+0.05
  d = det(G_379)
  if (d!=0) break 
  
}
print (i) ####### 1
print (d)## 3.424355e-53

########################################################

########################################################################################
Ginv <- solve(G_379)
Ginv[1:5,1:5]
dim(Ginv) # 567 567
###################################################################
source('C:/Users/uqsyada1/OneDrive - The University of Queensland/PhD/Trainings and workshop/201910_AsREML/Asreml-R workshop/Asreml(4-5)/Material_UQ45/full2sparse.R')
Ginv.sparse<-full2sparse(Ginv)
colnames(Ginv.sparse)<-c('Row','Column','Ainverse')
head(Ginv.sparse)

##################################################

class(Z_379_cl) ## character
cl_379 <- as.data.frame(Z_379_cl)
colnames(cl_379)[1]<-"clone"



######################################################
ahatinv <- Ginv.sparse
ahatinv<-as.matrix.data.frame(ahatinv)
attr(ahatinv,"rowNames")<-as.character(cl_379$clone)
attr(ahatinv,"colNames")<-as.character(cl_379$clone)
attr(ahatinv,"INVERSE")<-TRUE
head(ahatinv)

######################################################
library(asreml)
##################### TCH ######################################################


modelGBLUP_A = asreml(fixed=TCHBlup ~ Series + Region + Trial + Crop,
                      random=~vm(Clone,ahatinv),
                      #random=~Clone,
                      workspace=128e06,na.action = na.method(y="include"),
                      data= train_data)
summary(modelGBLUP_A)$varcomp
# #
# component std.error  z.ratio bound %ch
# vm(Clone, ahatinv)  92.96465  8.204566 11.33084     P   0
# units!R             42.06056  1.060154 39.67401     P   0

source('C:/Users/uqsyada1/OneDrive - The University of Queensland/PhD/Trainings and workshop/201910_AsREML/Asreml-R workshop/ASReml(4-5)/Material_UQ45/pin.r')
pin(modelGBLUP_A,h2_1~V1/(V1+V2))
# #
# Estimate         SE
# h2_1 0.6884984 0.01983252
##########################################################################

##############################################
BLUPS_TCH <- summary(modelGBLUP_A,coef=TRUE)$coef.random
BLUPS_TCH = as.data.frame(BLUPS_TCH)
head(BLUPS_TCH)
dim(BLUPS_TCH) # 379    3
a.effects <- BLUPS_TCH[1:411,1]
#######################################################

# getting prediction accuracies, i.e. correlating actual with predicted performance

all.names = c()
for (i in 1:nrow (BLUPS_TCH) ) {
  
  name = strsplit (rownames(BLUPS_TCH[i,]), split = "_") [[1]][2]
  
  all.names = c(all.names, name)
  
}

BLUPS_TCH$clone = all.names
head(BLUPS_TCH)

test_set$GEBV_TCH = 0

for (i in 1:nrow(test_set)){
  
  clone = noquote(test_set[i,'Clone'])
  
  sum <-BLUPS_TCH$solution[BLUPS_TCH$clone == clone][1] 
  
  test_set$GEBV_TCH[test_set$Clone == clone] = sum
  
}
###################################################################################
head(test_set,5)
cor (test_set$TCHBlup, test_set$GEBV_TCH,use = "pairwise.complete.obs")## 0.14,0.18
plot(test_set$TCHBlup, test_set$GEBV_TCH)

############################################################
## CCS

exclude <- intersect(train_data$Clone, test_set$Clone) ## 113
train_data$CCSBlup [train_data$Clone %in% exclude] = NA


modelGBLUP_A = asreml(fixed=CCSBlup ~ Series + Region + Trial + Crop,
                      random=~vm(Clone,ahatinv),
                      #random=~Clone,
                      workspace=128e06,na.action = na.method(y="include"),
                      data= train_data)
summary(modelGBLUP_A)$varcomp
# #
# component   std.error  z.ratio bound %ch
# vm(Clone, ahatinv) 0.4128673 0.034695726 11.89966     P   0
# units!R            0.0802726 0.002021174 39.71582     P   0

source('C:/Users/uqsyada1/OneDrive - The University of Queensland/PhD/Trainings and workshop/201910_AsREML/Asreml-R workshop/ASReml(4-5)/Material_UQ45/pin.r')
pin(modelGBLUP_A,h2_1~V1/(V1+V2))

# ##
# Estimate         SE
# h2_1 0.8372214 0.01199151
##########################################################################

##############################################
BLUPS_CCS <- summary(modelGBLUP_A,coef=TRUE)$coef.random
BLUPS_CCS = as.data.frame(BLUPS_CCS)
head(BLUPS_CCS)
dim(BLUPS_CCS) # 379    3
a.effects <- BLUPS_CCS[1:379,1]
#######################################################

# getting prediction accuracies, i.e. correlating actual with predicted performance

all.names = c()
for (i in 1:nrow (BLUPS_CCS) ) {
  
  name = strsplit (rownames(BLUPS_CCS[i,]), split = "_") [[1]][2]
  
  all.names = c(all.names, name)
  
}

BLUPS_CCS$clone = all.names
head(BLUPS_CCS)

test_set$GEBV_CCS = 0

for (i in 1:nrow(test_set)){
  
  clone = noquote(test_set[i,'Clone'])
  
  sum <-BLUPS_CCS$solution[BLUPS_CCS$clone == clone][1] 
  
  test_set$GEBV_CCS[test_set$Clone == clone] = sum
  
}
###################################################################################
head(test_set,5)
cor (test_set$CCSBlup, test_set$GEBV_CCS)## 0.2233292, 0.4386351
plot(test_set$CCSBlup, test_set$GEBV_CCS)

############################################################
### Fibre 

exclude <- intersect(train_data$Clone, test_set$Clone) ## 113
train_data$FibreBlup [train_data$Clone %in% exclude] = NA


modelGBLUP_A = asreml(fixed=FibreBlup ~ Series + Region + Trial + Crop,
                      random=~vm(Clone,ahatinv),
                      #random=~Clone,
                      workspace=128e06,na.action = na.method(y="include"),
                      data= train_data)
summary(modelGBLUP_A)$varcomp

# ## 
# component   std.error  z.ratio bound %ch
# vm(Clone, ahatinv)  2.447453 0.202202004 12.10400     P   0
# units!R             0.157607 0.003970085 39.69864     P   0

source('C:/Users/uqsyada1/OneDrive - The University of Queensland/PhD/Trainings and workshop/201910_AsREML/Asreml-R workshop/ASReml(4-5)/Material_UQ45/pin.r')
pin(modelGBLUP_A,h2_1~V1/(V1+V2))

# ##
# Estimate          SE
# h2_1 0.9394997 0.004919247
##########################################################################

##############################################
BLUPS_Fibre <- summary(modelGBLUP_A,coef=TRUE)$coef.random
BLUPS_Fibre = as.data.frame(BLUPS_Fibre)
head(BLUPS_Fibre)
dim(BLUPS_Fibre) # 567    3
a.effects <- BLUPS_Fibre[1:379,1]
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
cor (test_set$FibreBlup, test_set$GEBV_Fibre,use = "pairwise.complete.obs")## 0.3989093, 0.5459054
plot(test_set$FibreBlup, test_set$GEBV_Fibre)

################################################################
############################################################
##########################################################
### Region C
# Subset the data for the current region
region_data <- subset(Data, Region == "C")
region_data <- droplevels(region_data)
str(region_data) ## 201 LEVELS
# Set a seed for reproducibility
set.seed(42)

# Get the unique clones
unique_clones <- unique(region_data$Clone)
# Generate random indices for the test set
test_indices <- sample(length(unique_clones), floor(0.2 * length(unique_clones)))
# Create the test set by subsetting the data based on the selected clones and setting the response variable to NA
test_set <- region_data[region_data$Clone %in% unique_clones[test_indices], ]
test_set <- droplevels(test_set)
str(test_set) ## 40 clones

train_data <- region_data
exclude <- intersect(train_data$Clone, test_set$Clone) ## 82
train_data$TCHBlup [train_data$Clone %in% exclude] = NA

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
##################### TCH ######################################################


modelGBLUP_A = asreml(fixed=TCHBlup ~ Series + Region + Trial + Crop,
                      random=~vm(Clone,ahatinv),
                      #random=~Clone,
                      workspace=128e06,na.action = na.method(y="include"),
                      data= train_data)
summary(modelGBLUP_A)$varcomp
source('C:/Users/uqsyada1/OneDrive - The University of Queensland/PhD/Trainings and workshop/201910_AsREML/Asreml-R workshop/ASReml(4-5)/Material_UQ45/pin.r')
pin(modelGBLUP_A,h2_1~V1/(V1+V2))
##########################################################################

##############################################
BLUPS_TCH <- summary(modelGBLUP_A,coef=TRUE)$coef.random
BLUPS_TCH = as.data.frame(BLUPS_TCH)
head(BLUPS_TCH)
dim(BLUPS_TCH) # 411    3
a.effects <- BLUPS_TCH[1:201,1]
#######################################################

# getting prediction accuracies, i.e. correlating actual with predicted performance

all.names = c()
for (i in 1:nrow (BLUPS_TCH) ) {
  
  name = strsplit (rownames(BLUPS_TCH[i,]), split = "_") [[1]][2]
  
  all.names = c(all.names, name)
  
}

BLUPS_TCH$clone = all.names
head(BLUPS_TCH)

test_set$GEBV_TCH = 0

for (i in 1:nrow(test_set)){
  
  clone = noquote(test_set[i,'Clone'])
  
  sum <-BLUPS_TCH$solution[BLUPS_TCH$clone == clone][1] 
  
  test_set$GEBV_TCH[test_set$Clone == clone] = sum
  
}
###################################################################################
head(test_set,5)
cor (test_set$TCHBlup, test_set$GEBV_TCH,use = "pairwise.complete.obs")## 0.1053125, 0.1381576, -0.01463427
plot(test_set$TCHBlup, test_set$GEBV_TCH)

############################################################
## CCS

exclude <- intersect(train_data$Clone, test_set$Clone) ## 113
train_data$CCSBlup [train_data$Clone %in% exclude] = NA


modelGBLUP_A = asreml(fixed=CCSBlup ~ Series + Region + Trial + Crop,
                      random=~vm(Clone,ahatinv),
                      #random=~Clone,
                      workspace=128e06,na.action = na.method(y="include"),
                      data= train_data)
summary(modelGBLUP_A)$varcomp
source('C:/Users/uqsyada1/OneDrive - The University of Queensland/PhD/Trainings and workshop/201910_AsREML/Asreml-R workshop/ASReml(4-5)/Material_UQ45/pin.r')
pin(modelGBLUP_A,h2_1~V1/(V1+V2))
##########################################################################

##############################################
BLUPS_CCS <- summary(modelGBLUP_A,coef=TRUE)$coef.random
BLUPS_CCS = as.data.frame(BLUPS_CCS)
head(BLUPS_CCS)
dim(BLUPS_CCS) # 201    3
a.effects <- BLUPS_CCS[1:201,1]
#######################################################

# getting prediction accuracies, i.e. correlating actual with predicted performance

all.names = c()
for (i in 1:nrow (BLUPS_CCS) ) {
  
  name = strsplit (rownames(BLUPS_CCS[i,]), split = "_") [[1]][2]
  
  all.names = c(all.names, name)
  
}

BLUPS_CCS$clone = all.names
head(BLUPS_CCS)

test_set$GEBV_CCS = 0

for (i in 1:nrow(test_set)){
  
  clone = noquote(test_set[i,'Clone'])
  
  sum <-BLUPS_CCS$solution[BLUPS_CCS$clone == clone][1] 
  
  test_set$GEBV_CCS[test_set$Clone == clone] = sum
  
}
###################################################################################
head(test_set,5)
cor (test_set$CCSBlup, test_set$GEBV_CCS)## 0.4093551, 0.2217085, 0.3273817
plot(test_set$CCSBlup, test_set$GEBV_CCS)

############################################################
### Fibre 

exclude <- intersect(train_data$Clone, test_set$Clone) ## 113
train_data$FibreBlup [train_data$Clone %in% exclude] = NA


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
cor (test_set$FibreBlup, test_set$GEBV_Fibre,use = "pairwise.complete.obs")## 0.3235673, 0.1469458, 0.2054209
plot(test_set$FibreBlup, test_set$GEBV_Fibre)




