
######### Linkage-Disequilibrium Algorithm

## First calculate LD estimates among all pairs of SNPs using SelectionTools"

## Exactly like File1_LD_decay

library(data.table)
## reading the ld estimates

Z <- as.data.frame(fread("ld.txt",header=TRUE)) ## ld.txt is in karen new genetic map folder
LD <- Z
class(LD) ##  "data.frame"
rm(Z)
dim(LD) #### 340226655  6 (we have 26086 markers, pair wise r2 estimates)

#####################################################################################################
LD <- LD[,4:6]
head(LD)
max= max(LD$LD) ## 0.974117
min=min(LD$LD)  ## 0
################################################################################

nms = unique(c(LD$Name1,LD$Name2))
length(nms) # 26086

## convert table into matrix, where assigning 0 to pair-wise estimates to itself

mt = matrix(0,nrow=length(nms),ncol=length(nms))
dim(mt) ## 26086 26086

row.names(mt) = nms
colnames(mt) = nms

n=length(nms) ### 26086

offset=0

## setting a matrix

for (i in 1:26085){ #n
  #i=1 
  mt[i,(i+1):n] = LD$LD[(offset+1):(offset+ (n-i))]
  offset=offset+(n-i)
}

############################################

## LD matrix is ready; BUT it is upper triangular
## use linear algebra to make its full symmetric matrix

mt_sym <- mt + t(mt) ## symmetric matrix

mt_sym[1:10,1:10]
######### read the genetic map
## for sugarcane only 4502 snps were mapped on the genome

## we already filterd the common snps from the SNP array and snps which were positioned on the genome

df1 <- as.data.frame(fread("fil_gmap_4502.txt",header=TRUE))
dim(df1)

##############################################################
## this step is to find out how many common snps we have

markername <- colnames(mt_sym)
commonn <- subset(markername, markername%in%df1$Marker)

length(commonn) # 4502
head(df1)
####################################################################

#######################################################
## The genetic map with only common SNPs

df_m <- df1[df1$Marker%in% commonn,]
df_m <- droplevels(df_m)
dim(df_m) ## 4502    3
head(df_m) 
################
## 20-fold cross validation
## preparing the data in a matrix form
## rows represent the unmapped markers (for cross-validation, masked mapped snps)
## column represent the mapped markers

data1 <- mt_sym[(markername %in% commonn),(markername %in% commonn)]
class(data1)

dim(data1) ## 4502 4502 (for sugarcane)

################################################################
r = rep(seq(1,20), 225) ######### preparing the data for 20 fold cross-validation
r = c(r,1:2)
length(r)
table(r)
df_m$no = r
head(df_m)

###################### storing results ##################################################
res1 = c()
res2 = c()
res3 = c()
res4 = c()
res5 = c()

table(df_m$no)

for( r in 1:20){
  df2 = df_m$Marker[df_m$no == r]
  df_m_ss = df_m[df_m$Marker %in% df2,]
  
  data2 <- data1[(rownames(data1) %in% df2), !(rownames(data1) %in% df2)]
  dim(data2) ## rows represents "unmapped" and column represents "mapped" markers
  
################# LD algorithm ########################
  marker.names.all = c()
  chrom = c()
  Pos = c()
  ###########@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
  ###########@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  for ( i in 1:nrow(data2)){
    data3 <- data2[i,] ## data2 is numeric
    
    data4 <- sort(data3, decreasing = TRUE)
    head(data4)
 ## Step 1: check after sorting, if first two ld estimates is greater than 0.1   
    if(sum(data4[1:2] > 0.1) == 2) ## check alternate threshold, depending on species
    { 
      marker.names=names(data4[1:2]) ## if yes, checking the names of these two snps
      checkpos1 = subset(df1, Marker== marker.names[1])## position of first snp
      checkpos2 = subset(df1, Marker == marker.names[2])## postion of second snp
   ## Step 2: checking if both snps are on the same chromosome
      
      if (checkpos1$Chromosome == checkpos2$Chromosome)
      {name = rownames(data2)[i] ## if yes, checking the name of unmapped snp
      marker.names.all = c(marker.names.all, name)
      ch <- checkpos1$Chromosome
      chrom= c(chrom,ch)
      ## first LD
      LD1 <- data1[name,marker.names[1]]
      LD2 <- data1[name,marker.names[2]]
      pos1 <- checkpos1$Position
      pos2 <- checkpos2$Position
  ## Step 3 : Assigning position
      
      Mpos <- ((pos1 *LD1) + (pos2 * LD2))/(LD1 + LD2) ## assigning position of unmapped snps where weights are assigned based on LD estimates
      Pos <- c(Pos, Mpos)
 }
    }
    
  }
  ## storing all the resutls
  
  res1_tmp = length(marker.names.all)
  res1_tmp = round(res1_tmp/nrow(df_m_ss) * 100 , 1)
  res1 = c(res1, res1_tmp) ### representing efficiency, how many markers we can mapped
  
  df_new <- data.frame(Marker = marker.names.all,Chromosome=chrom, Position=Pos)
  df_m_merge = merge(df_new, df_m_ss, by = "Marker") ## merging new dataframe with old dataframe (having marker name, chromosome and position)
  
  # write.csv(df_m_merge, paste0("results",r,".csv")) # write resutls for each cross validation
  res2_tmp = sum(df_m_merge$Chromosome.x == df_m_merge$Chromosome.y)
  res2_tmp = round (res2_tmp / nrow(df_m_merge) * 100 , 1)
  res2 = c(res2, res2_tmp) ## represents the accuracy, how many markers are on the right chromosome
  
  df_m_merge$chrom_check = df_m_merge$Chromosome.x == df_m_merge$Chromosome.y
  df_m_merge_ss = df_m_merge [df_m_merge$chrom_check == TRUE,]
  dim(df_m_merge_ss)
  
  res3_tmp = round(mean(abs(df_m_merge_ss$Position.y - df_m_merge_ss$Position.x)),2)
  res3 = c(res3, res3_tmp) ## mean of absolute difference of actual and new assigned position
  
  res4_tmp = round(median(abs(df_m_merge_ss$Position.y - df_m_merge_ss$Position.x)),2)
  res4 = c(res4, res4_tmp) ## median of absolute difference of actual and new assigned position
  
  res5_tmp = round(sd(abs(df_m_merge_ss$Position.y - df_m_merge_ss$Position.x)),2)
  res5 = c(res5, res5_tmp) ## standard deviation of absolute difference of actual and new assigned position
  
  
  print(paste0("done with rep ", r))
  
}
  