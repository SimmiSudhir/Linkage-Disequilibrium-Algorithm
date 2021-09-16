###########################################

########### LD decay #######################################

rm(list=ls())
library(data.table)
marker <- as.data.frame(fread("Sugarcane_MarkerData.txt",header=TRUE))

gmap <- as.data.frame(fread("fil_gmap_4502.txt",header=TRUE))

########## Markers which are mapped on the Genetic (Q208) genetic map
common <- which(marker$name %in% gmap$Marker)
marker <- marker[common,]
fwrite(marker, "mar_4502.txt",row.names=F,col.names=T, sep= " ")
#########################################

######## Intra-chromosomal LD estimates

###############

library(SelectionTools)
st.input.dir <- "input"
st.output.dir <- "output"
##########################################################################
########## Reading marker data and Genetic map #################################

st.read.marker.data("mar_test.txt", format= "m")
st.read.map("map_test.txt", m.stretch = 1, format = "mcp", skip = 1)

st.copy.marker.data("q3", source.data.set="default" )

## calculate LD for each pair of SNPs
ld <- st.calc.ld(ld.measure="r2", data.set="q3")
class(ld) ## data.frame
dim(ld) 
####################
#######

dfr <- read.table("map_test.txt", header=TRUE)
str(dfr)
dfr$Chromosome <- as.factor(dfr$Chromosome)
levels(dfr$Chromosome)

###############################################################################################
## set up an empty results table
## calculate distance for every pair of SNPs

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


res$LD <- ld$LD
head(res)

max(res$Distance) ## make it subset of distance having 100cM ## use geom_smooth (method="loess")
head(res)

##############################################################
####### plotting decay of LD with distance
## and add a smooth line to show the trend of decay
## codes below represents plot a moving average of LD in increments of 5 cM distances.
#####################
res <- subset(res, Distance <= 100)

bins= seq(1,max(res$Distance),by=0.01)

my.means=rep(0, length(bins))
LD.averages =data.frame(bins, my.means)
for( i in 1:length(bins)){
  data_interval=subset(res,(res$Distance >= bins[i] & res$Distance < (bins[i] + 0.5)))
  LD.averages$my.means[i]=mean(data_interval$LD)
}

head(LD.averages)
LD.averages$Crop <- "Sugarcane"
colnames(LD.averages)[1] <- "bins"
colnames(LD.averages)[2] <- "means"

library(ggplot2)
p1 <- ggplot(data = LD.averages, aes(x = bins, y = means)) +
  geom_point(alpha = 0.1,color="blue")+ geom_smooth(method="loess")+
  geom_hline(yintercept=0.2, linetype="dashed", color = "red")+
  geom_hline(yintercept=0.1, linetype="dashed", color = "red")+
  labs(title = "Decay of Linkage Disequilibrium",
       x = "Distance between SNPs (cM)",
       y = expression (~r^2)) +   
  theme_bw()+ ylim(0, 0.5)
############ 

################## another way of representation 

library(dplyr)
library(stringr)

res1 <- res
res2 <- subset(res1, Distance <= 100)

res2$distc <- cut(res2$Distance, breaks=seq(from=min(res2$Distance)-1,to=max(res2$Distance)+1,by=20))

levels(res2$distc) <- c(levels(res2$distc),"None")
res2$distc[is.na(res2$distc)] <- "None"
res2$distc <- gsub("None", "> 99",res2$distc)

res3 <- res2 %>% group_by(distc) %>% summarise(mean=mean(LD),median=median(LD))
res3 <- res3 %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                        mid=start+((end-start)/2))

ggplot()+
  geom_point(data=res3,aes(x=start,y=mean),size=0.4,colour="grey20")+
  geom_line(data=res3,aes(x=start,y=mean),size=0.3,alpha=0.5,colour="grey40")+
  labs(x="Distance (cM)",y=expression(LD~(r^{2})))+
  theme_bw()

#####################################################

res <- subset(res, Distance <= 100)

bins= seq(1,max(res$Distance),by=0.01)

my.means=rep(0, length(bins))
LD.averages =data.frame(bins, my.means)
for( i in 1:length(bins)){
  data_interval=subset(res,(res$Distance >= bins[i] & res$Distance < (bins[i] + 0.5)))
  LD.averages$my.means[i]=mean(data_interval$LD)
}

head(LD.averages)
LD.averages$Crop <- "Sugarcane"
colnames(LD.averages)[1] <- "bins"
colnames(LD.averages)[2] <- "means"

library(ggplot2)
p1 <- ggplot(data = LD.averages, aes(x = bins, y = means)) +
  geom_point(alpha = 0.1,color="blue")+ geom_smooth(method="loess")+
  geom_hline(yintercept=0.2, linetype="dashed", color = "red")+
  geom_hline(yintercept=0.1, linetype="dashed", color = "red")+
  labs(title = "Decay of Linkage Disequilibrium",
       x = "Distance between SNPs (cM)",
       y = expression (~r^2)) +   
  theme_bw()+ ylim(0, 0.5)

###################### similarly combining all data together on one graph

data1 <- as.data.frame(fread("LD_averages_Sugarcane.txt",header=TRUE))
data2 <- as.data.frame(fread("LD_averages_wheat.txt",header=TRUE))
data3 <- as.data.frame(fread("LD_averages_Barley.txt",header=TRUE))
data4 <- rbind(data1,data2,data3)

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

ggsave("ld_decaytiff", units="in", width=6, height=3, dpi=300, compression = 'lzw')

dev.off()





