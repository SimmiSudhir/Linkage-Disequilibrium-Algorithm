rm(list=ls())
library(data.table)
data <- as.data.frame(fread("pheno_1318.txt", header=TRUE))
str(data)
data$Series <- as.factor(data$Series)
data$Region <- as.factor(data$Region)
data$Trial <- as.factor(data$Trial)
data$Crop <- as.factor(data$Crop)
data$Clone <- as.factor(data$Clone)
str(data)


# ########################################################
# 
# Series	A	  C	  N	   S	 Total
# 2013	  4	  1	  2	   3	 10
# 2014	 284	46	134	145	 609
# 2015	 256	91	122	119	 588
# 2016	  3	  62	130	167	362
# 2017	 53	  6	  21	36	116

########################################################
### exclude 2013 year data set and Region "C"
Data <- subset(data, Series != 2013 & Region != "C")

Data <- droplevels(Data)
str(Data) ## 1200 levels of clones with 4 levels of Series


clones_1200 <- unique(Data$Clone)
clones_1200 <- as.data.frame(clones_1200)
colnames(clones_1200)[1] <- "name"
#fwrite(clones_1200, "clones_1200.txt",col.names=T,row.names=F,sep = " ")

###########################################
## Read z2909 matrix

Z <- as.data.frame(fread("Marker_1318_dip.txt",header=T))
dim(Z) ## 1318 25753
Z[1:5,1:5]

clone_1318 <- as.data.frame(fread("clone_1318.txt", header=TRUE))
head(clone_1318)
Z_1318 <- cbind(clone_1318,Z)
Z_1318[1:5,1:5]

Z_1200 <- Z_1318[Z_1318$x %in% clones_1200$name,]
Z_1200[1:5,1:5]

Z1200 <- Z_1200[-1]
row.names(Z1200) <- Z_1200$x
snp_matrix <- as.matrix(Z1200)
########################################
# ## make GRM 
# library(AGHmatrix)
# G <- Gmatrix(snp_matrix, missingValue=-9, 
#              maf=0.01, method="VanRaden",ploidy=2)
## make GRM 
library(BGLR)
X <- scale(snp_matrix, center=TRUE, scale=TRUE)
dim(X) ## 1200 25753
p <- ncol(X)
## COMPUTING g
G <- tcrossprod(X)/p
G[1:5,1:5]

dim(G) ## 1200 BY 1200
G[1:5,1:5]
GA <- G
dim(G) ## 1200 BY 1200
G[1:5,1:5]
GA <- G
det(GA) ## 0
#########################################################################
for (i in 1:100) {  # i added 8x 0.05 in total
  
  diag(GA)=diag(GA)+0.05
  d = det(GA)
  if (d!=0) break 
  
}
print (i) ####### 1
print (d) ##  3.077377e-318

########################################################################################
Ginv <- solve(GA)
Ginv[1:5,1:5]
dim(Ginv) # 1200 1200
###################################################################
source('C:/Users/uqsyada1/OneDrive - The University of Queensland/PhD/Trainings and workshop/201910_AsREML/Asreml-R workshop/Asreml(4-5)/Material_UQ45/full2sparse.R')
Ginv.sparse<-full2sparse(Ginv)
colnames(Ginv.sparse)<-c('Row','Column','Ainverse')
head(Ginv.sparse)

#fwrite (Ginv.sparse,"GRM_1200_dip.txt",col.names=T,row.names=F, sep=" ")

#################################################

ahatinv<-Ginv.sparse

ahatinv<-as.matrix.data.frame(ahatinv)
attr(ahatinv,"rowNames")<-as.character(Z_1200$x)
attr(ahatinv,"colNames")<-as.character(Z_1200$x)
attr(ahatinv,"INVERSE")<-TRUE
head(ahatinv)

#####################################################
Pheno <- Data

Year = 2017
dim(Pheno) # 35955     8
testpop = Pheno[Pheno$Series==Year,]
testpop <- droplevels(testpop)
nrow(testpop) #  4529
head(testpop)
dim(testpop) # 4529    8
str(testpop) # 691 levels of clones
str(Pheno) ## 2909 levels of clones
###################################################################################
trainpop <- Pheno
exclude = noquote(intersect ( trainpop$Clone [trainpop$Series == Year],  trainpop$Clone [trainpop$Series != Year]))
length(exclude) # 179 common clones; should be removed from the training population

trainpop$TCHBlup [trainpop$Series %in% Year] = NA
trainpop$TCHBlup [trainpop$Clone %in% exclude] = NA

trainpop$CCSBlup [trainpop$Series %in% Year] = NA
trainpop$CCSBlup [trainpop$Clone %in% exclude] = NA

trainpop$FibreBlup [trainpop$Series %in% Year] = NA
trainpop$FibreBlup [trainpop$Clone %in% exclude] = NA

######################################

head(trainpop)
# # Get the row names of snp_matrix
# snp_row_names <- rownames(snp_matrix)
# 
# # Reorder the 'Clone' column in trainpop based on the row names
# trainpop_reordered <- trainpop[match(snp_row_names, trainpop$Clone), ]


library(asreml)

modelGBLUP_A = asreml(fixed=TCHBlup ~ Series + Crop + Trial + Region,
                      random=~vm(Clone,ahatinv),
                      #random=~Clone,
                      workspace=128e06,na.action = na.method(y="include"),
                      data= trainpop)
summary(modelGBLUP_A)$varcomp

source('C:/Users/uqsyada1/OneDrive - The University of Queensland/PhD/Trainings and workshop/201910_AsREML/Asreml-R workshop/ASReml(4-5)/Material_UQ45/pin.r')
pin(modelGBLUP_A,h2_1~V1/(V1+V2))
##########################################################################
BLUPS_TCH <- summary(modelGBLUP_A,coef=TRUE)$coef.random
BLUPS_TCH = as.data.frame(BLUPS_TCH)
head(BLUPS_TCH)
#############################################################################
######################################################################
dim(BLUPS_TCH) #  1317    3
a.effects <- BLUPS_TCH[1:1200,1]
##########################################################################
# getting prediction accuracies, i.e. correlating actual with predicted performance

all.names = c()
for (i in 1:nrow (BLUPS_TCH) ) {
  
  name = strsplit (rownames(BLUPS_TCH[i,]), split = "_") [[1]][2]
  
  all.names = c(all.names, name)
  
}

BLUPS_TCH$clone = all.names
head(BLUPS_TCH)

testpop$GEBV_TCH = 0

for (i in 1:nrow(testpop)){
  
  clone = noquote(testpop[i,'Clone'])
  
  sum <-BLUPS_TCH$solution[BLUPS_TCH$clone == clone][1] 
  
  testpop$GEBV_TCH[testpop$Clone == clone] = sum
  
}
###################################################################################
head(testpop,5)
############### test accrding to region wise ##########################

cor (testpop$TCHBlup, testpop$GEBV_TCH) ##0.3789938


####################################################
########## CCS ###############################
library(asreml)

modelGBLUP_A = asreml(fixed=CCSBlup ~ Series+Crop+Trial + Region,
                      random=~vm(Clone,ahatinv),
                      #random=~Clone,
                      workspace=128e06,na.action = na.method(y="include"),
                      data= trainpop)
summary(modelGBLUP_A)$varcomp

source('C:/Users/uqsyada1/OneDrive - The University of Queensland/PhD/Trainings and workshop/201910_AsREML/Asreml-R workshop/ASReml(4-5)/Material_UQ45/pin.r')
pin(modelGBLUP_A,h2_1~V1/(V1+V2))

BLUPS_CCS <- summary(modelGBLUP_A,coef=TRUE)$coef.random
BLUPS_CCS = as.data.frame(BLUPS_CCS)
head(BLUPS_CCS)
#############################################################################
######################################################################
dim(BLUPS_CCS) #  1317    3
a.effects <- BLUPS_CCS[1:1200,1]

##########################################################################
# getting prediction accuracies, i.e. correlating actual with predicted performance

all.names = c()
for (i in 1:nrow (BLUPS_CCS) ) {
  
  name = strsplit (rownames(BLUPS_CCS[i,]), split = "_") [[1]][2]
  
  all.names = c(all.names, name)
  
}

BLUPS_CCS$clone = all.names
head(BLUPS_CCS)

testpop$GEBV_CCS = 0

for (i in 1:nrow(testpop)){
  
  clone = noquote(testpop[i,'Clone'])
  
  sum <-BLUPS_CCS$solution[BLUPS_CCS$clone == clone][1] 
  
  testpop$GEBV_CCS[testpop$Clone == clone] = sum
  
}
####################
cor (testpop$CCSBlup, testpop$GEBV_CCS) ## 0.2063132

######################################################

########## Fibre ###############################
library(asreml)

modelGBLUP_A = asreml(fixed=FibreBlup ~ Series+Crop+Trial + Region,
                      random=~vm(Clone,ahatinv),
                      #random=~Clone,
                      workspace=128e06,na.action = na.method(y="include"),
                      data= trainpop)
summary(modelGBLUP_A)$varcomp

source('C:/Users/uqsyada1/OneDrive - The University of Queensland/PhD/Trainings and workshop/201910_AsREML/Asreml-R workshop/ASReml(4-5)/Material_UQ45/pin.r')
pin(modelGBLUP_A,h2_1~V1/(V1+V2))

BLUPS_Fibre <- summary(modelGBLUP_A,coef=TRUE)$coef.random
BLUPS_Fibre = as.data.frame(BLUPS_Fibre)
head(BLUPS_Fibre)
#############################################################################
######################################################################
dim(BLUPS_Fibre) #  1317    3
a.effects <- BLUPS_Fibre[1:1200,1]

##########################################################################
# getting prediction accuracies, i.e. correlating actual with predicted performance

all.names = c()
for (i in 1:nrow (BLUPS_Fibre) ) {
  
  name = strsplit (rownames(BLUPS_Fibre[i,]), split = "_") [[1]][2]
  
  all.names = c(all.names, name)
  
}

BLUPS_Fibre$clone = all.names
head(BLUPS_Fibre)

testpop$GEBV_Fibre = 0

for (i in 1:nrow(testpop)){
  
  clone = noquote(testpop[i,'Clone'])
  
  sum <-BLUPS_Fibre$solution[BLUPS_Fibre$clone == clone][1] 
  
  testpop$GEBV_Fibre[testpop$Clone == clone] = sum
  
}
####################
cor (testpop$FibreBlup, testpop$GEBV_Fibre) ## 0.5648064


######################################################
#########################################################
## 26K continuous 
##################################################
rm(list=ls())
library(data.table)
data <- as.data.frame(fread("pheno_1318.txt", header=TRUE))
str(data)
data$Series <- as.factor(data$Series)
data$Region <- as.factor(data$Region)
data$Trial <- as.factor(data$Trial)
data$Crop <- as.factor(data$Crop)
data$Clone <- as.factor(data$Clone)
str(data)


# ########################################################
# 
# Series	A	  C	  N	   S	 Total
# 2013	  4	  1	  2	   3	 10
# 2014	 284	46	134	145	 609
# 2015	 256	91	122	119	 588
# 2016	  3	  62	130	167	362
# 2017	 53	  6	  21	36	116

########################################################
### exclude 2013 year data set and Region "C"
Data <- subset(data, Series != 2013 & Region != "C")

Data <- droplevels(Data)
str(Data) ## 1200 levels of clones with 4 levels of Series


clones_1200 <- unique(Data$Clone)
clones_1200 <- as.data.frame(clones_1200)
colnames(clones_1200)[1] <- "name"
#fwrite(clones_1200, "clones_1200.txt",col.names=T,row.names=F,sep = " ")

###########################################
## Read z2909 matrix

Z <- as.data.frame(fread("Marker_1318_con.txt",header=T))
dim(Z) ## 1318 25753
Z[1:5,1:5]

clone_1318 <- as.data.frame(fread("clone_1318.txt", header=TRUE))
head(clone_1318)
Z_1318 <- cbind(clone_1318,Z)
Z_1318[1:5,1:5]

Z_1200 <- Z_1318[Z_1318$x %in% clones_1200$name,]
Z_1200[1:5,1:5]

Z1200 <- Z_1200[-1]
row.names(Z1200) <- Z_1200$x
snp_matrix <- as.matrix(Z1200)
########################################
## make GRM 
library(BGLR)
X <- scale(snp_matrix, center=TRUE, scale=TRUE)
dim(X) ## 1200 25753
p <- ncol(X)
## COMPUTING g
G <- tcrossprod(X)/p
G[1:5,1:5]

dim(G) ## 1200 BY 1200
G[1:5,1:5]
GA <- G
det(GA) ## 0
#########################################################################
for (i in 1:100) {  # i added 8x 0.05 in total
  
  diag(GA)=diag(GA)+0.05
  d = det(GA)
  if (d!=0) break 
  
}
print (i) ####### 2
print (d) ##  1.276848e-306

########################################################################################
Ginv <- solve(GA)
Ginv[1:5,1:5]
dim(Ginv) # 1200 1200
###################################################################
source('C:/Users/uqsyada1/OneDrive - The University of Queensland/PhD/Trainings and workshop/201910_AsREML/Asreml-R workshop/Asreml(4-5)/Material_UQ45/full2sparse.R')
Ginv.sparse<-full2sparse(Ginv)
colnames(Ginv.sparse)<-c('Row','Column','Ainverse')
head(Ginv.sparse)

#fwrite (Ginv.sparse,"GRM_1200_dip.txt",col.names=T,row.names=F, sep=" ")

#################################################

ahatinv<-Ginv.sparse

ahatinv<-as.matrix.data.frame(ahatinv)
attr(ahatinv,"rowNames")<-as.character(Z_1200$x)
attr(ahatinv,"colNames")<-as.character(Z_1200$x)
attr(ahatinv,"INVERSE")<-TRUE
head(ahatinv)

#####################################################
Pheno <- Data

Year = 2017
dim(Pheno) # 15157     8
testpop = Pheno[Pheno$Series==Year,]
testpop <- droplevels(testpop)
nrow(testpop) # 684
head(testpop)
dim(testpop) # 684   8
str(testpop) # 94 levels of clones
str(Pheno) ## 1200 levels of clones
###################################################################################
trainpop <- Pheno
exclude = noquote(intersect ( trainpop$Clone [trainpop$Series == Year],  trainpop$Clone [trainpop$Series != Year]))
length(exclude) # 179 common clones; should be removed from the training population

trainpop$TCHBlup [trainpop$Series %in% Year] = NA
trainpop$TCHBlup [trainpop$Clone %in% exclude] = NA

trainpop$CCSBlup [trainpop$Series %in% Year] = NA
trainpop$CCSBlup [trainpop$Clone %in% exclude] = NA

trainpop$FibreBlup [trainpop$Series %in% Year] = NA
trainpop$FibreBlup [trainpop$Clone %in% exclude] = NA

######################################

head(trainpop)
# # Get the row names of snp_matrix
# snp_row_names <- rownames(snp_matrix)
# 
# # Reorder the 'Clone' column in trainpop based on the row names
# trainpop_reordered <- trainpop[match(snp_row_names, trainpop$Clone), ]


library(asreml)

modelGBLUP_A = asreml(fixed=TCHBlup ~ Series + Crop + Trial + Region,
                      random=~vm(Clone,ahatinv),
                      #random=~Clone,
                      workspace=128e06,na.action = na.method(y="include"),
                      data= trainpop)
summary(modelGBLUP_A)$varcomp

source('C:/Users/uqsyada1/OneDrive - The University of Queensland/PhD/Trainings and workshop/201910_AsREML/Asreml-R workshop/ASReml(4-5)/Material_UQ45/pin.r')
pin(modelGBLUP_A,h2_1~V1/(V1+V2))
##########################################################################
BLUPS_TCH <- summary(modelGBLUP_A,coef=TRUE)$coef.random
BLUPS_TCH = as.data.frame(BLUPS_TCH)
head(BLUPS_TCH)
#############################################################################
######################################################################
dim(BLUPS_TCH) #  1200    3
a.effects <- BLUPS_TCH[1:1200,1]
##########################################################################
# getting prediction accuracies, i.e. correlating actual with predicted performance

all.names = c()
for (i in 1:nrow (BLUPS_TCH) ) {
  
  name = strsplit (rownames(BLUPS_TCH[i,]), split = "_") [[1]][2]
  
  all.names = c(all.names, name)
  
}

BLUPS_TCH$clone = all.names
head(BLUPS_TCH)

testpop$GEBV_TCH = 0

for (i in 1:nrow(testpop)){
  
  clone = noquote(testpop[i,'Clone'])
  
  sum <-BLUPS_TCH$solution[BLUPS_TCH$clone == clone][1] 
  
  testpop$GEBV_TCH[testpop$Clone == clone] = sum
  
}
###################################################################################
head(testpop,5)
############### test accrding to region wise ##########################

cor (testpop$TCHBlup, testpop$GEBV_TCH) ##0.3855197


####################################################
########## CCS ###############################
library(asreml)

modelGBLUP_A = asreml(fixed=CCSBlup ~ Series+Crop+Trial + Region,
                      random=~vm(Clone,ahatinv),
                      #random=~Clone,
                      workspace=128e06,na.action = na.method(y="include"),
                      data= trainpop)
summary(modelGBLUP_A)$varcomp

source('C:/Users/uqsyada1/OneDrive - The University of Queensland/PhD/Trainings and workshop/201910_AsREML/Asreml-R workshop/ASReml(4-5)/Material_UQ45/pin.r')
pin(modelGBLUP_A,h2_1~V1/(V1+V2))

BLUPS_CCS <- summary(modelGBLUP_A,coef=TRUE)$coef.random
BLUPS_CCS = as.data.frame(BLUPS_CCS)
head(BLUPS_CCS)
#############################################################################
######################################################################
dim(BLUPS_CCS) #  1200    3
a.effects <- BLUPS_CCS[1:1200,1]

##########################################################################
# getting prediction accuracies, i.e. correlating actual with predicted performance

all.names = c()
for (i in 1:nrow (BLUPS_CCS) ) {
  
  name = strsplit (rownames(BLUPS_CCS[i,]), split = "_") [[1]][2]
  
  all.names = c(all.names, name)
  
}

BLUPS_CCS$clone = all.names
head(BLUPS_CCS)

testpop$GEBV_CCS = 0

for (i in 1:nrow(testpop)){
  
  clone = noquote(testpop[i,'Clone'])
  
  sum <-BLUPS_CCS$solution[BLUPS_CCS$clone == clone][1] 
  
  testpop$GEBV_CCS[testpop$Clone == clone] = sum
  
}
####################
cor (testpop$CCSBlup, testpop$GEBV_CCS) ## 0.1817881

 ##########################
####################################

########## Fibre ###############################
library(asreml)

modelGBLUP_A = asreml(fixed=FibreBlup ~ Series+Crop+Trial + Region,
                      random=~vm(Clone,ahatinv),
                      #random=~Clone,
                      workspace=128e06,na.action = na.method(y="include"),
                      data= trainpop)
summary(modelGBLUP_A)$varcomp

source('C:/Users/uqsyada1/OneDrive - The University of Queensland/PhD/Trainings and workshop/201910_AsREML/Asreml-R workshop/ASReml(4-5)/Material_UQ45/pin.r')
pin(modelGBLUP_A,h2_1~V1/(V1+V2))

BLUPS_Fibre <- summary(modelGBLUP_A,coef=TRUE)$coef.random
BLUPS_Fibre = as.data.frame(BLUPS_Fibre)
head(BLUPS_Fibre)
#############################################################################
######################################################################
dim(BLUPS_Fibre) #  1200    3
a.effects <- BLUPS_Fibre[1:1200,1]

##########################################################################
# getting prediction accuracies, i.e. correlating actual with predicted performance

all.names = c()
for (i in 1:nrow (BLUPS_Fibre) ) {
  
  name = strsplit (rownames(BLUPS_Fibre[i,]), split = "_") [[1]][2]
  
  all.names = c(all.names, name)
  
}

BLUPS_Fibre$clone = all.names
head(BLUPS_Fibre)

testpop$GEBV_Fibre = 0

for (i in 1:nrow(testpop)){
  
  clone = noquote(testpop[i,'Clone'])
  
  sum <-BLUPS_Fibre$solution[BLUPS_Fibre$clone == clone][1] 
  
  testpop$GEBV_Fibre[testpop$Clone == clone] = sum
  
}
####################
cor (testpop$FibreBlup, testpop$GEBV_Fibre) ## 0.5644059



###########################################################
########################################################
## continuyous 58k

##################################################
rm(list=ls())
library(data.table)
data <- as.data.frame(fread("pheno_1318.txt", header=TRUE))
str(data)
data$Series <- as.factor(data$Series)
data$Region <- as.factor(data$Region)
data$Trial <- as.factor(data$Trial)
data$Crop <- as.factor(data$Crop)
data$Clone <- as.factor(data$Clone)
str(data)


# ########################################################
# 
# Series	A	  C	  N	   S	 Total
# 2013	  4	  1	  2	   3	 10
# 2014	 284	46	134	145	 609
# 2015	 256	91	122	119	 588
# 2016	  3	  62	130	167	362
# 2017	 53	  6	  21	36	116

########################################################
### exclude 2013 year data set and Region "C"
Data <- subset(data, Series != 2013 & Region != "C")

Data <- droplevels(Data)
str(Data) ## 1200 levels of clones with 4 levels of Series


clones_1200 <- unique(Data$Clone)
clones_1200 <- as.data.frame(clones_1200)
colnames(clones_1200)[1] <- "name"
#fwrite(clones_1200, "clones_1200.txt",col.names=T,row.names=F,sep = " ")

###########################################
## Read z2909 matrix

Z <- as.data.frame(fread("Marker_1318_con_58364.txt",header=T))
dim(Z) ## 1318 25753
Z[1:5,1:5]

clone_1318 <- as.data.frame(fread("clone_1318.txt", header=TRUE))
head(clone_1318)
#Z_1318 <- cbind(clone_1318,Z)
#Z_1318[1:5,1:5]

Z_1200 <- Z[Z$Clone %in% clones_1200$name,]
Z_1200[1:5,1:5]

Z1200 <- Z_1200[-1]
row.names(Z1200) <- Z_1200$Clone
snp_matrix <- as.matrix(Z1200)
########################################
## make GRM 
library(BGLR)
X <- scale(snp_matrix, center=TRUE, scale=TRUE)
dim(X) ## 1200 25753
p <- ncol(X)
## COMPUTING g
G <- tcrossprod(X)/p
G[1:5,1:5]

dim(G) ## 1200 BY 1200
G[1:5,1:5]
GA <- G
det(GA) ## 0
#########################################################################
for (i in 1:100) {  # i added 8x 0.05 in total
  
  diag(GA)=diag(GA)+0.05
  d = det(GA)
  if (d!=0) break 
  
}
print (i) ####### 2
print (d) ##  3.200029e-290

########################################################################################
Ginv <- solve(GA)
Ginv[1:5,1:5]
dim(Ginv) # 1200 1200
###################################################################
source('C:/Users/uqsyada1/OneDrive - The University of Queensland/PhD/Trainings and workshop/201910_AsREML/Asreml-R workshop/Asreml(4-5)/Material_UQ45/full2sparse.R')
Ginv.sparse<-full2sparse(Ginv)
colnames(Ginv.sparse)<-c('Row','Column','Ainverse')
head(Ginv.sparse)

#fwrite (Ginv.sparse,"GRM_1200_dip.txt",col.names=T,row.names=F, sep=" ")

#################################################

ahatinv<-Ginv.sparse

ahatinv<-as.matrix.data.frame(ahatinv)
attr(ahatinv,"rowNames")<-as.character(Z_1200$Clone)
attr(ahatinv,"colNames")<-as.character(Z_1200$Clone)
attr(ahatinv,"INVERSE")<-TRUE
head(ahatinv)

#####################################################
Pheno <- Data

Year = 2017
dim(Pheno) # 15157     8
testpop = Pheno[Pheno$Series==Year,]
testpop <- droplevels(testpop)
nrow(testpop) # 684
head(testpop)
dim(testpop) # 684   8
str(testpop) # 94 levels of clones
str(Pheno) ## 1200 levels of clones
###################################################################################
trainpop <- Pheno
exclude = noquote(intersect ( trainpop$Clone [trainpop$Series == Year],  trainpop$Clone [trainpop$Series != Year]))
length(exclude) # 179 common clones; should be removed from the training population

trainpop$TCHBlup [trainpop$Series %in% Year] = NA
trainpop$TCHBlup [trainpop$Clone %in% exclude] = NA

trainpop$CCSBlup [trainpop$Series %in% Year] = NA
trainpop$CCSBlup [trainpop$Clone %in% exclude] = NA

trainpop$FibreBlup [trainpop$Series %in% Year] = NA
trainpop$FibreBlup [trainpop$Clone %in% exclude] = NA

######################################
head(trainpop)

library(asreml)

modelGBLUP_A = asreml(fixed=TCHBlup ~ Series + Crop + Trial + Region,
                      random=~vm(Clone,ahatinv),
                      #random=~Clone,
                      workspace=128e06,na.action = na.method(y="include"),
                      data= trainpop)
summary(modelGBLUP_A)$varcomp

source('C:/Users/uqsyada1/OneDrive - The University of Queensland/PhD/Trainings and workshop/201910_AsREML/Asreml-R workshop/ASReml(4-5)/Material_UQ45/pin.r')
pin(modelGBLUP_A,h2_1~V1/(V1+V2))
##########################################################################
BLUPS_TCH <- summary(modelGBLUP_A,coef=TRUE)$coef.random
BLUPS_TCH = as.data.frame(BLUPS_TCH)
head(BLUPS_TCH)
#############################################################################
######################################################################
dim(BLUPS_TCH) #  1200    3
a.effects <- BLUPS_TCH[1:1200,1]
##########################################################################
# getting prediction accuracies, i.e. correlating actual with predicted performance

all.names = c()
for (i in 1:nrow (BLUPS_TCH) ) {
  
  name = strsplit (rownames(BLUPS_TCH[i,]), split = "_") [[1]][2]
  
  all.names = c(all.names, name)
  
}

BLUPS_TCH$clone = all.names
head(BLUPS_TCH)

testpop$GEBV_TCH = 0

for (i in 1:nrow(testpop)){
  
  clone = noquote(testpop[i,'Clone'])
  
  sum <-BLUPS_TCH$solution[BLUPS_TCH$clone == clone][1] 
  
  testpop$GEBV_TCH[testpop$Clone == clone] = sum
  
}
###################################################################################
head(testpop,5)
############### test accrding to region wise ##########################

cor (testpop$TCHBlup, testpop$GEBV_TCH) ##0.3954051

####################################################
########## CCS ###############################
library(asreml)

modelGBLUP_A = asreml(fixed=CCSBlup ~ Series+Crop+Trial + Region,
                      random=~vm(Clone,ahatinv),
                      #random=~Clone,
                      workspace=128e06,na.action = na.method(y="include"),
                      data= trainpop)
summary(modelGBLUP_A)$varcomp

source('C:/Users/uqsyada1/OneDrive - The University of Queensland/PhD/Trainings and workshop/201910_AsREML/Asreml-R workshop/ASReml(4-5)/Material_UQ45/pin.r')
pin(modelGBLUP_A,h2_1~V1/(V1+V2))

BLUPS_CCS <- summary(modelGBLUP_A,coef=TRUE)$coef.random
BLUPS_CCS = as.data.frame(BLUPS_CCS)
head(BLUPS_CCS)
#############################################################################
######################################################################
dim(BLUPS_CCS) #  1200    3
a.effects <- BLUPS_CCS[1:1200,1]

##########################################################################
# getting prediction accuracies, i.e. correlating actual with predicted performance

all.names = c()
for (i in 1:nrow (BLUPS_CCS) ) {
  
  name = strsplit (rownames(BLUPS_CCS[i,]), split = "_") [[1]][2]
  
  all.names = c(all.names, name)
  
}

BLUPS_CCS$clone = all.names
head(BLUPS_CCS)

testpop$GEBV_CCS = 0

for (i in 1:nrow(testpop)){
  
  clone = noquote(testpop[i,'Clone'])
  
  sum <-BLUPS_CCS$solution[BLUPS_CCS$clone == clone][1] 
  
  testpop$GEBV_CCS[testpop$Clone == clone] = sum
  
}
####################
cor (testpop$CCSBlup, testpop$GEBV_CCS) ## 0.1953652

######################################################

########## Fibre ###############################
library(asreml)

modelGBLUP_A = asreml(fixed=FibreBlup ~ Series+Crop+Trial + Region,
                      random=~vm(Clone,ahatinv),
                      #random=~Clone,
                      workspace=128e06,na.action = na.method(y="include"),
                      data= trainpop)
summary(modelGBLUP_A)$varcomp

source('C:/Users/uqsyada1/OneDrive - The University of Queensland/PhD/Trainings and workshop/201910_AsREML/Asreml-R workshop/ASReml(4-5)/Material_UQ45/pin.r')
pin(modelGBLUP_A,h2_1~V1/(V1+V2))

BLUPS_Fibre <- summary(modelGBLUP_A,coef=TRUE)$coef.random
BLUPS_Fibre = as.data.frame(BLUPS_Fibre)
head(BLUPS_Fibre)
#############################################################################
######################################################################
dim(BLUPS_Fibre) #  1200    3
a.effects <- BLUPS_Fibre[1:1200,1]

##########################################################################
# getting prediction accuracies, i.e. correlating actual with predicted performance

all.names = c()
for (i in 1:nrow (BLUPS_Fibre) ) {
  
  name = strsplit (rownames(BLUPS_Fibre[i,]), split = "_") [[1]][2]
  
  all.names = c(all.names, name)
  
}

BLUPS_Fibre$clone = all.names
head(BLUPS_Fibre)

testpop$GEBV_Fibre = 0

for (i in 1:nrow(testpop)){
  
  clone = noquote(testpop[i,'Clone'])
  
  sum <-BLUPS_Fibre$solution[BLUPS_Fibre$clone == clone][1] 
  
  testpop$GEBV_Fibre[testpop$Clone == clone] = sum
  
}
####################
cor (testpop$FibreBlup, testpop$GEBV_Fibre) ##0.5620017


