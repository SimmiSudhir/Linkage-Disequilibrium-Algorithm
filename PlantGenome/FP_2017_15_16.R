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
Data <- subset(data, Series != 2013 & Series != 2014 & Region != "C")

Data <- droplevels(Data)
str(Data) ## 771 levels of clones with 4 levels of Series

clones_771 <- unique(Data$Clone)
clones_771 <- as.data.frame(clones_771)
colnames(clones_771)[1] <- "name"

Z <- as.data.frame(fread("Marker_1318_dip.txt",header=T))
dim(Z) ## 1318 25753
Z[1:5,1:5]

clone_1318 <- as.data.frame(fread("clone_1318.txt", header=TRUE))
head(clone_1318)
Z_1318 <- cbind(clone_1318,Z)
Z_1318[1:5,1:5]

Z_771 <- Z_1318[Z_1318$x %in% clones_771$name,]
Z_771[1:5,1:5]

Z771 <- Z_771[-1]
row.names(Z771) <- Z_771$x
snp_matrix <- as.matrix(Z771)
## make GRM 
library(AGHmatrix)
G <- Gmatrix(snp_matrix, missingValue=-9, 
             maf=0.01, method="VanRaden",ploidy=2)

dim(G) ## 771 771
G[1:5,1:5]
GA <- G
det(GA) ## -2.682177e-247

for (i in 1:100) {  # i added 8x 0.05 in total
  
  diag(GA)=diag(GA)+0.05
  d = det(GA)
  if (d!=0) break 
  
}
print (i) ####### 1
print (d) # 3.378056e-188

#######################################################################################
Ginv <- solve(GA)
Ginv[1:5,1:5]
dim(Ginv) # 1200 1200
###################################################################
source('C:/Users/uqsyada1/OneDrive - The University of Queensland/PhD/Trainings and workshop/201910_AsREML/Asreml-R workshop/Asreml(4-5)/Material_UQ45/full2sparse.R')
Ginv.sparse<-full2sparse(Ginv)
colnames(Ginv.sparse)<-c('Row','Column','Ainverse')
head(Ginv.sparse)

fwrite (Ginv.sparse,"GRM_771_dip.txt",col.names=T,row.names=F, sep=" ")

###################################################

ahatinv<-read.table("GRM_771_dip.txt",h=T)

ahatinv<-as.matrix.data.frame(ahatinv)
attr(ahatinv,"rowNames")<-as.character(Z_771$x)
attr(ahatinv,"colNames")<-as.character(Z_771$x)
attr(ahatinv,"INVERSE")<-TRUE
head(ahatinv)
###############################################
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

#############################################
trainpop <- Pheno
exclude = noquote(intersect ( trainpop$Clone [trainpop$Series == Year],  trainpop$Clone [trainpop$Series != Year]))
length(exclude) # 179 common clones; should be removed from the training population

trainpop$TCHBlup [trainpop$Series %in% Year] = NA
trainpop$TCHBlup [trainpop$Clone %in% exclude] = NA

trainpop$CCSBlup [trainpop$Series %in% Year] = NA
trainpop$CCSBlup [trainpop$Clone %in% exclude] = NA

trainpop$FibreBlup [trainpop$Series %in% Year] = NA
trainpop$FibreBlup [trainpop$Clone %in% exclude] = NA
##########################################################
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
dim(BLUPS_TCH) #  1317    3
a.effects <- BLUPS_TCH[1:771,1]
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

cor (testpop$TCHBlup, testpop$GEBV_TCH) ##0.428

############### test accrding to region wise ##########################
test = testpop[testpop$Region == "N",]
cor (test$TCHBlup, test$GEBV_TCH) #0.0408772

#########################################################################
################################################################################
test = testpop[testpop$Region == "A",]
cor (test$TCHBlup, test$GEBV_TCH) ##0.5091918

#############################################

####################################################
test = testpop[testpop$Region == "S",]
cor (test$TCHBlup, test$GEBV_TCH)  ## 0.36
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
a.effects <- BLUPS_CCS[1:771,1]

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
cor (testpop$CCSBlup, testpop$GEBV_CCS) ## 0.2936104

############### test accrding to region wise ##########################
test = testpop[testpop$Region == "N",]
cor (test$CCSBlup, test$GEBV_CCS) # 0.3026186

#########################################################################
################################################################################
test = testpop[testpop$Region == "A",]
cor (test$CCSBlup, test$GEBV_CCS) ## 0.2453181

#############################################

####################################################
test = testpop[testpop$Region == "S",]
cor (test$CCSBlup, test$GEBV_CCS)  ## 0.4061933

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
a.effects <- BLUPS_Fibre[1:771,1]

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
cor (testpop$FibreBlup, testpop$GEBV_Fibre) ## 0.5286548

############### test accrding to region wise ##########################
test = testpop[testpop$Region == "N",]
cor (test$FibreBlup, test$GEBV_Fibre) # 0.603445

#########################################################################
################################################################################
test = testpop[testpop$Region == "A",]
cor (test$FibreBlup, test$GEBV_Fibre) ## 0.4513917

#############################################

####################################################
test = testpop[testpop$Region == "S",]
cor (test$FibreBlup, test$GEBV_Fibre)  ## 0.6337471

######################################################


