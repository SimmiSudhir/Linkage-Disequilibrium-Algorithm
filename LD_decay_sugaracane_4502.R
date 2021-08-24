rm(list=ls())
library(data.table)
marker <- as.data.frame(fread("genos_SelTools.txt",header=TRUE))
gmap <- as.data.frame(fread("fil_gmap_4502.txt",header=TRUE))

common <- which(marker$name %in% gmap$Marker)
marker <- marker[common,]
fwrite(marker, "mar_4502.txt",row.names=F,col.names=T, sep= " ")
# #cc <- intersect(marker$name,gmap$Marker)
# length(cc)
# marker1 <- marker[marker$name %in% cc,]
# dim



####################################################################################
library(SelectionTools)
st.input.dir <- "input"
st.output.dir <- "output"
##########################################################################
########## Reading marker data and Genetic map #################################

st.read.marker.data("mar_4502.txt", format= "m")
st.read.map("fil_gmap_4502.txt", m.stretch = 1, format = "mcp", skip = 1)## 

st.copy.marker.data("q3", source.data.set="default" )

## calculate LD
ld <- st.calc.ld(ld.measure="r2", data.set="q3")
class(ld) ## data.frame
dim(ld) ## 209635      6
 # length(levels(as.factor(ld$Chrom))) ## 129
nms = unique(c(ld$Name1,ld$Name2))
# #length(nms) ## 4373 (## make sense 4502-4373= 129 chromosome)
# head(ld)

#########################################################################
str(ld)
ld$Chrom <- as.factor(ld$Chrom)
length(levels(ld$Chrom)) ## 129
################################################################################################
dfr <- read.table("fil_gmap_4502.txt", header=TRUE)
str(dfr)
dfr$Chromosome <- as.factor(dfr$Chromosome)
levels(dfr$Chromosome)
# dfr_49 <- subset(dfr, Chromosome==49) ## maximum marker at this chromosome

###############################################################################################
res= c()
for(r in levels(dfr$Chromosome))
{
  dfr_ss = subset(dfr, Chromosome == r)
  dfr_ss <- droplevels(dfr_ss)
  
  my.results <- data.frame(matrix(0, nrow=choose(nrow(dfr_ss),2),ncol=1))
  colnames(my.results) <- "Distance"
  counter = 1
  for ( i in 1: (nrow(dfr_ss)-1)){
    
    for (j in (i+1):nrow(dfr_ss)) {
      
      my.results$Distance[counter]= abs(dfr_ss$Position[j]-dfr_ss$Position[i])
      
      counter=counter+1
      
    }
  }
  res = rbind(res, my.results)
  
}

## dim(res) ## 209635      1
## dim(ld) ## 209635      6

res$LD <- ld$LD
head(res)
dim(res) ## 209635     2
# View(res) 
max(res$Distance) ## make it subset of distance having 50cm ## use geom_smooth (method="loess")
head(res)

# ########################################################################
### practice 
res1 <- res
res2 <- subset(res1, Distance <= 100)

res2$distc <- cut(res2$Distance, breaks=seq(from=min(res2$Distance)-1,to=max(res2$Distance)+1,by=20))
# library(tidyr) (tidyr is not installing)
levels(res2$distc) <- c(levels(res2$distc),"None")
res2$distc[is.na(res2$distc)] <- "None"
res2$distc <- gsub("None", "> 99",res2$distc)
library(dplyr)
library(stringr)
res3 <- res2 %>% group_by(distc) %>% summarise(mean=mean(LD),median=median(LD))
res3 <- res3 %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                        mid=start+((end-start)/2))

ggplot()+
  geom_point(data=res3,aes(x=start,y=mean),size=0.4,colour="grey20")+
  geom_line(data=res3,aes(x=start,y=mean),size=0.3,alpha=0.5,colour="grey40")+
  labs(x="Distance (Megabases)",y=expression(LD~(r^{2})))+
  theme_bw()

library(ggplot2)
p1<- ggplot(res2,aes(x=distc,y=LD))
p1 <- p1+geom_boxplot()
p1
##

res <- subset(res, Distance <= 100) ## 209635 obs. of  2 variables:
# dim(res) ## 72355     2
# max(res$Distance) ## 50
# 
# ##########################################################################
library(ggplot2)
p <- ggplot(data = res, aes(x = Distance, y = LD)) +
  geom_point(alpha = 0.1,color="blue")+
  labs(title = "Decay of Linkage Disequilibrium",
       x = "Distance between SNPs (cM)",
       y = "expression (r^2)") +
  theme_bw()

################################################################
bins= seq(1,max(res$Distance),by=0.01)

my.means=rep(0, length(bins))
LD.averages =data.frame(bins, my.means)
for( i in 1:length(bins)){
  data_interval=subset(res,(res$Distance >= bins[i] & res$Distance < (bins[i] + 0.01)))
  LD.averages$my.means[i]=mean(data_interval$LD)
}

# ################################################################
# p <- p+geom_point(data=LD.averages,aes(x=LD.averages$bins,y=LD.averages$my.means),col="red")
# p <- p + geom_hline(yintercept=0.1, linetype="dashed", color = "black")
# #################################################################
head(LD.averages)
LD.averages$Crop <- "Sugarcane"
colnames(LD.averages)[1] <- "bins"
colnames(LD.averages)[2] <- "means"
fwrite(LD.averages,file='LD_averages_Sugarcane.txt',row.names=F, col.names=T,sep= " ")
##################################################
###################################################################


library(ggplot2)
p1 <- ggplot(data = LD.averages, aes(x = bins, y = means)) +
  geom_point(alpha = 0.1,color="blue")+ geom_smooth(method="loess")+
   geom_hline(yintercept=0.2, linetype="dashed", color = "red")+
  geom_hline(yintercept=0.1, linetype="dashed", color = "red")+
  labs(title = "Decay of Linkage Disequilibrium",
       x = "Distance between SNPs (cM)",
       y = expression (~r^2)) +   
  theme_bw()+ ylim(0, 0.5)

###################################################################
##################wheat #####################################

library(SelectionTools)
st.input.dir <- "input" # Input files
st.output.dir <- "output" # Output files
st.read.marker.data ("wheat_genos.txt" ,format="m") ## 460 individuals with 37861 markers
st.restrict.marker.data ( NoAll.MAX=2 )  ## 460 ind and 37861 markers
st.restrict.marker.data ( MaMis.MAX=0.1 ) ##460 ind and 28794 markers

st.restrict.marker.data ( ExHet.MIN= 0.1 ) ## 460 ind and 18475 markers 
st.restrict.marker.data ( InMis.MAX=0.1 ) ## 450  ind and 18475 markers
st.read.map ("wheat_gmap.txt", format="mcp", skip=1)

st.copy.marker.data ( "q3", "default" )
ld_w <- st.calc.ld ( ld.measure="r2", 
                   data.set="q3" )

dim(ld_w) ## 10985846        6

#########################################################################
str(ld_w)
ld_w$Chrom <- as.factor(ld_w$Chrom)
length(levels(ld_w$Chrom)) ## 21
################################################################################################


###################################################################
X = st.marker.data.statistics(data.set="default")
Marker_names <- data.frame (Marker = X$marker.list$Name)

##########################################3333333
# w_marker <- as.data.frame(fread("wheat_genos.txt",header=TRUE))
w_gmap <- as.data.frame(fread("wheat_gmap.txt",header=TRUE))
com <- intersect(w_gmap$name , Marker_names$Marker)
length(com) ## 18475
w_gmap <- w_gmap[w_gmap$name%in%com,] ## genetic map of only 18475 markers
head(w_gmap)
dim(w_gmap)
str(w_gmap) ## 18475     3
w_gmap$chr <- as.factor(w_gmap$chr)
levels(w_gmap$chr)  ## 21

table(w_gmap$chr)
#######################################################
w_gmap$pos <- w_gmap$pos/100
head(w_gmap)
## genetic map markers are not in the same order of ld_w
## make it in the same order
od_gmap <- w_gmap[match(Marker_names$Marker,w_gmap$name),]
head(od_gmap)
od_gmap$chr <- as.factor(od_gmap$chr)
########################################################################
w_res= c()
for(r in levels(od_gmap$chr))
{
  w_gmap_ss = subset(od_gmap, chr == r)
  w_gmap_ss <- droplevels(w_gmap_ss)
  
  w_my.results <- data.frame(matrix(0, nrow=choose(nrow(w_gmap_ss),2),ncol=1))
  colnames(w_my.results) <- "Distance"
  counter = 1
  for ( i in 1: (nrow(w_gmap_ss)-1)){
    
    for (j in (i+1):nrow(w_gmap_ss)) {
      
      w_my.results$Distance[counter]= abs(w_gmap_ss$pos[j]-w_gmap_ss$pos[i])
      
      counter=counter+1
      
    }
  }
  w_res = rbind(w_res, w_my.results)
  
} ## start at 5.14 pm on August 18, 2021 ##
############# it will take time simi
############# Be patient
############ Be patient
dim(w_res) ## 10985846 1
w_res$LD <- ld_w$LD
head(w_res) ## 10985846 1
wres<- w_res
### practice 
#res <- subset(res, Distance <= 100)
wres$distc <- cut(wres$Distance, breaks=seq(from=min(wres$Distance),to=max(wres$Distance),by=5))
library(ggplot2)
p1<- ggplot(wres,aes(x=distc,y=LD))
p1 <- p1+geom_boxplot()
p1







# w_res$Distance <- (w_res$Distance)/100
max(w_res$Distance) ## 234.7

w_Res <- subset(w_res, Distance <= 100) ## 10305924        2
fwrite(w_Res,file='wheat_dis_ld_18475.txt',row.names=F, col.names=T,sep= " ")
bins_w= seq(1,max(w_Res$Distance),by=0.5)

my.means_w=rep(0, length(bins_w))
LD.averages_w =data.frame(bins_w, my.means_w)
for( i in 1:length(bins_w)){
  data_interval_w=subset(w_Res,(w_Res$Distance >= bins_w[i] & w_Res$Distance < (bins_w[i] + 0.5)))
  LD.averages_w$my.means_w[i]=mean(data_interval_w$LD)
}

LD.averages_w$Crop <- "Wheat"
colnames(LD.averages_w)[1] <- "bins"
colnames(LD.averages_w)[2] <- "means"
fwrite(LD.averages_w,file='LD_averages_wheat.txt',row.names=F, col.names=T,sep= " ")
#############################################################################
library(ggplot2)

p2 <- ggplot(data = data3, aes(x = bins, y = means)) +
  geom_point(alpha = 0.1,color="blue")+ geom_smooth(method="loess")+
  geom_hline(yintercept=0.2, linetype="dashed", color = "red")+
  geom_hline(yintercept=0.1, linetype="dashed", color = "red")+
  labs(title = "Decay of Linkage Disequilibrium",
       x = "Distance between SNPs (cM)",
       y = expression (~r^2)) +   
  theme_bw()+ ylim(0, 0.5)
##############################################################################
################## Barley####################################################
####### need to include barley ld decay graph on the same graph ##############




################################################################
########## I am thinking to combine both data with different lines
rm(list=ls())
library(data.table)

data1 <- as.data.frame(fread("LD_averages_Sugarcane.txt",header=TRUE))
data2 <- as.data.frame(fread("LD_averages_wheat.txt",header=TRUE))
data3 <- as.data.frame(fread("LD_averages_Barley.txt",header=TRUE))
data4 <- rbind(data1,data2,data3)
# head(data4[data4$bins<3,])
#install.packages("ggpubr")
# library(ggpubr)
library(ggplot2)
p1 <- ggplot(data = data4, aes(x = bins, y = means,color=Crop)) +
  geom_point(alpha=0.3) +
  geom_smooth(aes( color=Crop),span=0.5)+
  scale_color_manual(values=c("red","blue","green"))+
  geom_hline(yintercept=0.1, linetype="dashed", color = "black")+
  labs(title = " LD Decay",
       x = "Distance between SNPs (cM)",
       y = expression (~r^2)) +   
  theme_bw()+ ylim(0, 0.5)+
  theme(legend.position = "right", legend.text = element_text(size=16))+
theme(axis.title.x = element_text( size=16)) +
  theme(axis.title.y = element_text( size=16))+
  theme(strip.text.x = element_text(size = 12))
  p1 + xlim(c(0,10))
 
ggsave("ld_decay.tiff", units="in", width=6, height=3, dpi=300, compression = 'lzw')

dev.off()
######################################################################################
install.packages("viridis")
library(viridis)
p2 <- ggplot(data = data4, aes(x = bins, y = means,color=Crop, fill=Crop)) +
  geom_point(alpha=0.3) +geom_smooth(aes( color=Crop),method="loess",se=FALSE)+
  scale_color_viridis(direction = 1,discrete=TRUE,option="C")+scale_fill_viridis(direction = 1,discrete=TRUE,option="C")+
  geom_hline(yintercept=0.1, linetype="dashed", color = "black")+
  labs(title = " LD Decay",
       x = "Distance between SNPs (cM)",
       y = expression (~r^2)) +   
  theme_bw()+ ylim(0, 0.5)+
  theme(legend.position = "right", legend.text = element_text(size=13))+
  theme(axis.title.x = element_text( size=13)) +
  theme(axis.title.y = element_text( size=13))+
  theme(strip.text.x = element_text(size = 12))

ggsave("ld_decay_withoutSE.tiff", units="in", width=6, height=3, dpi=300, compression = 'lzw')

dev.off()













