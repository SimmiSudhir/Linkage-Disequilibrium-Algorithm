##BayesR ## further Prediciton Scenario 1

## training population 2015,2015 and 2016 and predict 2017
## remove year 2013 and Region C
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


### exclude 2013 year data set and Region "C"
Data <- subset(data, Series != 2013 & Region != "C")

Data <- droplevels(Data)
str(Data) ## 1200 levels of clones with 4 levels of Series
## 3 levels of Region

########################################################
library(asreml)

model = asreml(fixed=TCHBlup ~ Series + Region + Trial + Crop +Clone,
               workspace=128e06,na.action = na.method(y="include"),
               data= Data)
summary(model)
Blues_tch <- summary(model,coef=TRUE)$coef.fixed
class(Blues_tch) ## data.frame
head(Blues_tch)
dim(Blues_tch) ## 1263    3

Blues_tch <- Blues_tch[1:1200,]
Blues_tch <- as.data.frame(Blues_tch)
dim(Blues_tch) ## 1200    3
head(Blues_tch)
tail(Blues_tch)
##########################
new_df <- data.frame(
  id = sub("^Clone_", "", row.names(Blues_tch)),
  solution = Blues_tch$solution
  
)
colnames(new_df)[2] <- "TCH_BLUEs"
###############################################
model1 = asreml(fixed=CCSBlup ~ Series + Region + Trial + Crop +Clone,
               workspace=128e06,na.action = na.method(y="include"),
               data= Data)
summary(model1)
Blues_CCS <- summary(model1,coef=TRUE)$coef.fixed
Blues_CCS <- Blues_CCS[1:1200,]
Blues_CCS <- as.data.frame(Blues_CCS)
head(Blues_CCS)

new_df <- cbind(new_df,Blues_CCS$solution)
colnames(new_df)[3] <-"CCS_BLUEs"

#####################################################
model2 = asreml(fixed=FibreBlup ~ Series + Region + Trial + Crop +Clone,
                workspace=128e06,na.action = na.method(y="include"),
                data= Data)
summary(model2)
Blues_Fibre <- summary(model2,coef=TRUE)$coef.fixed
Blues_Fibre <- Blues_Fibre[1:1200,]
Blues_Fibre <- as.data.frame(Blues_Fibre)
head(Blues_Fibre)

new_df <- cbind(new_df,Blues_Fibre$solution)
colnames(new_df)[4] <-"Fibre_BLUEs"
####################################################

head(new_df)
dim(new_df) ## 1200 4

fwrite(new_df,"BLUES_1200_Scenario1.txt", col.names=T, row.names=F, sep=" ")

#############################################################

##### make a training pheno where we need to predict 2017 clones
##################################################################
## make training and test population
############ all overlapping clones in training population
####### from test population should be NA ################

#####################################################
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
#######
###############################################################
######### Preparing training population ######################
#########
trainpop <- new_df
trainpop[trainpop$id %in% testpop$Clone, c(2,3,4)] <- NA
head(trainpop)

num_rows_with_na <- sum(apply(is.na(trainpop), 1, any)) ## 94

########################################################
#######preparing pheno file for gctb software
tch_pheno <- trainpop$TCH_BLUEs
tch_pheno <- as.data.frame(tch_pheno)
tch_pheno$gid <- seq(1:nrow(trainpop))
tch_pheno$famid <- seq(1:nrow(trainpop))

tch_pheno <- tch_pheno[,c(2,3,1)]
write.table(tch_pheno,"Scenario1_gt_tch.pheno.txt",col.names=F, row.names=F, sep=" ")

##########################################################
CCS_pheno <- trainpop$CCS_BLUEs
CCS_pheno <- as.data.frame(CCS_pheno)
CCS_pheno$gid <- seq(1:nrow(trainpop))
CCS_pheno$famid <- seq(1:nrow(trainpop))

CCS_pheno <- CCS_pheno[,c(2,3,1)]
write.table(CCS_pheno,"Scenario1_gt_ccs.pheno.txt",col.names=F, row.names=F, sep=" ")

#############################################################
Fibre_pheno <- trainpop$Fibre_BLUEs
Fibre_pheno <- as.data.frame(Fibre_pheno)
Fibre_pheno$gid <- seq(1:nrow(trainpop))
Fibre_pheno$famid <- seq(1:nrow(trainpop))

Fibre_pheno <- Fibre_pheno[,c(2,3,1)]
write.table(Fibre_pheno,"Scenario1_gt_fibre.pheno.txt",col.names=F, row.names=F, sep=" ")

##################################################################
########## genotypic data
library(data.table)
M_1318 <- as.data.frame(fread("Marker_1318_dip.txt", header=T))
M_1318[1:5,1:5]
dim(M_1318) # 1318 25753

M_1318_total <- as.data.frame(fread("Marker_1318_con_58364.txt", header=T))
dim(M_1318_total) ## 1318 58364
M_1318_total[1:5,1:5]

row.names(M_1318) <- M_1318_total$Clone

M_1318[1:5,1:5]

M_1200 <- M_1318[row.names(M_1318) %in% trainpop$id,]
dim(M_1200) ## 1200 25753
M_1200[1:5,1:5]


row_names_match <- all(rownames(M_1200) == trainpop$id)
new_M_1200 <- M_1200[match(trainpop$id, rownames(M_1200)), ]

Marker_dip_1200 <- new_M_1200
gt <- Marker_dip_1200
gt[gt == "0"] <- "A A"
gt[gt == "1"] <- "A T"
gt[gt == "2"] <- "T T"

gt[1:10,1:5]
class(gt)

new_df <- as.data.frame(matrix(0,nrow=nrow(gt),ncol=1))

names(new_df) <- "Fam_Id"
new_df$Fam_Id <- seq(1:nrow(new_df))
new_df$Ind_ID <- seq(1:nrow(new_df))
new_df$Pat_ID <- "0"
new_df$Mat_ID <- "0"
new_df$Sex <- "0"
new_df$Phenotype <- "-9"
new.gt <- cbind(new_df,gt)

########################################
# format: no headers

fwrite(new.gt, "Scenario1_gt.ped", row.names=F,col.names = F, sep=" ", quote=FALSE)
########
MNAMES <- colnames(M_1318)
head(MNAMES)

length(MNAMES) ## 25753
df <- as.data.frame(MNAMES)
df$Chromosome <- 1
df$Position <- seq(1:nrow(df))

gmap <- df[,c(2,1,3)]

write.table(gmap, "Scenario1_gt.map",
            col.names = FALSE, row.names = FALSE,
            quote = FALSE)
#######################################################