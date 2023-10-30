rm(list=ls())
setwd("C:/Users/uqsyada1/OneDrive - The University of Queensland/PhD/Research Chapter-Allele dosages/Data/Output/15Plates_data/total_batch_file/Analysis_Review/Final_Analysis/Folder1/BayesR_Scenario2/con_58k")
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
Data <- subset(data, Series != 2013 & Series != 2014 &Region != "C")

Data <- droplevels(Data)
str(Data)

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
dim(Pheno) #  10053     8
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
#################################################################

M_1318 <- as.data.frame(fread("Marker_1318_con_58364.txt",header=T))
dim(M_1318) ## 1318 58364
M_1318[1:5,1:5]
row.names(M_1318) <- M_1318$Clone
M_1318 <- M_1318[-1]

M_771 <- M_1318[row.names(M_1318) %in% trainpop$id,]
dim(M_771) ## 771 25753
M_771[1:5,1:5]

# Get the row names of M_771 and assign them to a variable
row_names <- row.names(M_771)

# Use the match() function to check the order
# Check if the order matches
order_match <- identical(row_names, trainpop$id)

# Print the order match
print(order_match) ## FALSE

# Check for mismatches
#mismatch_indices <- which(row_names != trainpop$id)

row_names_match <- all(rownames(M_771) == trainpop$id)
new_M_771 <- M_771[match(trainpop$id, rownames(M_771)), ]

Marker_con_771 <- new_M_771

dim(Marker_con_771) ## 1200 25753

Marker_con_771 <- as.matrix(Marker_con_771)
Marker_con_771 <- t(Marker_con_771)
dim(Marker_con_771)
Marker_con_771[1:5,1:5]

fwrite(Marker_con_771,"Scenario2.gt.cgen",row.names=F, col.names=F,sep=" ")

############ fam file ##########################

fam_df <- as.data.frame(matrix(0,nrow=nrow(M_771),ncol=1))


names(fam_df) <- "Fam_Id"
fam_df$Fam_Id <- seq(1:nrow(fam_df))
fam_df$Ind_ID <- seq(1:nrow(fam_df))
fam_df$Pat_ID <- "-9"
fam_df$Mat_ID <- "-9"
fam_df$Sex <- "-9"

fam_df$Pheno <- "-9"

head(fam_df)

fwrite(fam_df, "Scenario2.gt.fam", row.names=F,col.names = F, sep=" ", quote=FALSE)

################## bim file ##################################

MNAMES <- colnames(M_1318)
head(MNAMES)

length(MNAMES) ## 25753
df <- as.data.frame(MNAMES)
df$Chromosome <- 1
df$Position <- seq(1:nrow(df))
df$GD <- 0
head(df)
df <- df[,c(2,1,4,3)]
df1 <- df
df1$a1 <- "-9"
df1$a2 <- "-9"

head(df1)
fwrite(df1, "Scenario2.gt.bim",col.names = FALSE, row.names = FALSE,sep=" ",quote = FALSE)

#############################################################


