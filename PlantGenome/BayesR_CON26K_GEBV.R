###############################################
############ Scenario 2 Bayes R ############################
#############################################
###### 771 clones ##################################

###############################################
############ Scenario 1 Bayes R ############################
#############################################

###########First test population ##################
######################1200 clones####################
rm(list=ls())
setwd("C:/Users/uqsyada1/OneDrive - The University of Queensland/PhD/Research Chapter-Allele dosages/Data/Output/15Plates_data/total_batch_file/Analysis_Review/Final_Analysis/Folder1/BayesR_Scenario2/con_26k")
library(data.table)
data <- as.data.frame(fread("pheno_1318.txt", header=TRUE))
str(data)
data$Series <- as.factor(data$Series)
data$Region <- as.factor(data$Region)
data$Trial <- as.factor(data$Trial)
data$Crop <- as.factor(data$Crop)
data$Clone <- as.factor(data$Clone)
str(data)


### exclude 2013, 2014 year data set and Region "C"
Data <- subset(data, Series != 2013 & Series != 2014 & Region != "C")

Data <- droplevels(Data)
str(Data) ## 1200 levels of clones with 4 levels of Series
## 3 levels of Region

new_df <- as.data.frame(fread("BLUES_771_Scenario2.txt",header=TRUE))

Pheno <- Data

# Match and replace TCHBlup column
Pheno$TCHBlup <- new_df$TCH_BLUEs[match(Pheno$Clone, new_df$id)]

# Match and replace CCSBlup column
Pheno$CCSBlup <- new_df$CCS_BLUEs[match(Pheno$Clone, new_df$id)]

# Match and replace FibreBlup column
Pheno$FibreBlup <- new_df$Fibre_BLUEs[match(Pheno$Clone, new_df$id)]

head(Pheno)

Year = 2017
dim(Pheno) #  15157     8
testpop = Pheno[Pheno$Series==Year,]
testpop <- droplevels(testpop)
nrow(testpop) #  684
head(testpop)
dim(testpop) # 684   8
str(testpop) # 94 levels of clones
str(Pheno) ## 1200 levels of clones
################################################

##################################################
######### continuous 26k

############################################
## Read z matrix

M_1318 <- as.data.frame(fread("Marker_1318_con.txt",header=T))
dim(M_1318) ## 1318 25753
M_1318[1:5,1:5]

M_1318_total <- as.data.frame(fread("Marker_1318_con_58364.txt", header=T))
dim(M_1318_total) ## 1318 58364
M_1318_total[1:5,1:5]

row.names(M_1318) <- M_1318_total$Clone

M_771 <- M_1318[row.names(M_1318) %in% new_df$id,]
dim(M_771) ## 771 25753
M_771[1:5,1:5]

Z1200 <- M_771

snp_matrix <- as.matrix(Z1200)
#########################################################

## make GRM 
library(BGLR)
X <- scale(snp_matrix, center=TRUE, scale=TRUE)
#X <- snp_matrix
########## TCH ##############################

#snp_tch <- as.data.frame(fread("New_out.snpRes_26000_tch.txt",header=T))
snp_tch <- as.data.frame(fread("Scenario2_Out_TCH.snpRes.txt", header=T))

head(snp_tch)
tch_effects <- snp_tch$A1Effect
length(tch_effects) ## 25753
hist(tch_effects)

plot(tch_effects^2, ylab='Estimated Squared-Marker Effect',
     type='o',cex=.5,col=4,main='Marker Effects')

Effects <- as.matrix(tch_effects)
dim(Effects)
Effects[1:5,]

p_bv_tch <- X %*% Effects ## X matrix would be 1200 by 257
p_bv_tch<- as.data.frame(p_bv_tch)
colnames(p_bv_tch) <- "tch_gebv"
head(p_bv_tch)



testpop$tch_GEBV = 0

for( i in 1:nrow(testpop)){
  Clone = noquote(testpop[i, 'Clone'])
  sum <- p_bv_tch$tch_gebv[row.names(p_bv_tch) == Clone]
  testpop$tch_GEBV[testpop$Clone == Clone] =sum
  
}



cor(testpop$TCHBlup, testpop$tch_GEBV) ## 0.5140468

######################################################
######################## CCS #######################
snp_ccs <- as.data.frame(fread("Scenario2_Out_CCS.snpRes.txt", header=T))
head(snp_ccs)
ccs_effects <- snp_ccs$A1Effect
length(ccs_effects) ## 25753
hist(ccs_effects)

plot(ccs_effects^2, ylab='Estimated Squared-Marker Effect',
     type='o',cex=.5,col=4,main='Marker Effects')

CCS_Effects <- as.matrix(ccs_effects)
dim(CCS_Effects)
CCS_Effects[1:5,]

p_bv_ccs <- X %*% CCS_Effects ## X matrix would be 1200 by 257
p_bv_ccs<- as.data.frame(p_bv_ccs)
colnames(p_bv_ccs) <- "ccs_gebv"
head(p_bv_ccs)



testpop$ccs_GEBV = 0

for( i in 1:nrow(testpop)){
  Clone = noquote(testpop[i, 'Clone'])
  sum <- p_bv_ccs$ccs_gebv[row.names(p_bv_ccs) == Clone]
  testpop$ccs_GEBV[testpop$Clone == Clone] =sum
  
}



cor(testpop$CCSBlup, testpop$ccs_GEBV) ## 0.2135017
#######################################################
####################### Fibre ##########################
snp_fibre <- as.data.frame(fread("Scenario2_Out_Fibre.snpRes.txt", header=T))
head(snp_fibre)
fibre_effects <- snp_fibre$A1Effect
length(fibre_effects) ## 25753
hist(fibre_effects)

plot(fibre_effects^2, ylab='Estimated Squared-Marker Effect',
     type='o',cex=.5,col=4,main='Marker Effects')

fibre_Effects <- as.matrix(fibre_effects)
dim(fibre_Effects)
fibre_Effects[1:5,]

p_bv_fibre <- X %*% fibre_Effects ## X matrix would be 1200 by 257
p_bv_fibre<- as.data.frame(p_bv_fibre)
colnames(p_bv_fibre) <- "fibre_gebv"
head(p_bv_fibre)



testpop$fibre_GEBV = 0

for( i in 1:nrow(testpop)){
  Clone = noquote(testpop[i, 'Clone'])
  sum <- p_bv_fibre$fibre_gebv[row.names(p_bv_fibre) == Clone]
  testpop$fibre_GEBV[testpop$Clone == Clone] =sum
  
}

cor(testpop$FibreBlup, testpop$fibre_GEBV) ## 0.5854675
###############################################################
####################MSE ##################################
############ MSE #########################################

MSE_TCH <- sqrt(mean((testpop$TCHBlup - testpop$tch_GEBV)^2)) ## 
MSE_CCS <- sqrt(mean((testpop$CCSBlup - testpop$ccs_GEBV)^2)) ## 
MSE_Fibre <- sqrt(mean((testpop$FibreBlup - testpop$fibre_GEBV)^2))## 

####################################################
