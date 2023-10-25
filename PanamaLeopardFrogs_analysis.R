#load packages
library(vegan)
library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)
library(RColorBrewer)
###########################################################################
#############################Beta diversity################################
###########################################################################

#load in asv table
asv.tbl<-read.delim('PanamaLeopardFrogs/asv_table.txt', row.names = 1, header=T)

#load in meta data
meta<-read.delim('PanamaLeopardFrogs/filtered_map_Piedra.txt', header=T)

#subset asv table to remove non-Alo de Piedra samples (that's all that's in the filtered map)
asv.tbl<-asv.tbl[,names(asv.tbl) %in% meta$SampleID]

#get column sums for rarefaction
min(colSums(asv.tbl))
#1322

#rarefy table
asv.rare<-rrarefy(t(asv.tbl), sample=1322)

#calculate beta diversity
ko_pcoa<-capscale(asv.rare  ~ 1, distance='jaccard')

#pull out x/y coordinates
ko.scores<-scores(ko_pcoa)

#grab only sample coordinates, write to data frame
ko.coords<-as.data.frame(ko.scores$sites)

#create sample names as a column
ko.coords$SampleID<-row.names(ko.coords)

#map back meta data
ko.coords<-merge(ko.coords, meta, by.x='SampleID', by.y='SampleID')

#calculate percent variation explained for first two axis
100*round(ko_pcoa$CA$eig[1]/sum(ko_pcoa$CA$eig), 3)
#32.2
100*round(ko_pcoa$CA$eig[2]/sum(ko_pcoa$CA$eig), 3)
#12.6

#plot PCoA
library(ggplot2)
ggplot(ko.coords, aes(MDS1, MDS2, colour=Species))+
  geom_point(size=2)+
  scale_color_manual(values = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000'))+
  theme_bw()+
  guides(alpha = "none")+
  xlab("PC1- 32.2%")+
  ylab("PC2- 12.6%")
  
  
#calculate non-parametric permanova
#calculate distance
s16.dis<-as.data.frame(as.matrix(vegdist(asv.rare, method='jaccard')))
s16.dis2<-s16.dis
s16.dis2$SampleID<-row.names(s16.dis2)

#select only columns of interest for the metadata
meta2<-meta[,c(1,9)]

#add metadata 
s16.dis2<-merge(s16.dis2, meta2, by='SampleID')
  
#remove sample id column
s16.dis3<-s16.dis2[,2:91]
  
#calculate adonis
adonis2(s16.dis3[,-90] ~ Species, data=s16.dis3, permutations = 10000)

#          Df SumOfSqs      R2      F    Pr(>F)    
#Species   9  0.19297 0.38741 5.5511 9.999e-05 ***
#Residual 79  0.30515 0.61259                     
#Total    88  0.49812 1.00000         


###########################################################################
#############################Beta diversity################################
###########################################################################

#load in asv table
asv.tbl<-read.delim('PanamaLeopardFrogs/asv_table.txt', row.names = 1, header=T)

#subset asv table to remove non-Alo de Piedra samples (that's all that's in the filtered map)
asv.tbl<-asv.tbl[,names(asv.tbl) %in% meta$SampleID]

#load in meta data
meta<-read.delim('PanamaLeopardFrogs/filtered_map_Piedra.txt', header=T)

#CALCULATE RICHNESS & add metadata & statistics
larv.alph<-as.data.frame(specnumber(t(asv.tbl)))
larv.alph$SampleID<-row.names(larv.alph)
larv.alph<-merge(larv.alph, meta, by='SampleID')
larv.alph$Richness<-as.numeric(larv.alph$`specnumber(t(asv.tbl))`)
pairwise.t.test(larv.alph$Richness, larv.alph$Species, p.adjust.method = 'hochberg')
#no significant difference in richness, plot it anyways

#plot richness
ggplot(larv.alph, aes(Species, Richness, fill=Species))+
  geom_jitter()+
  geom_boxplot()+
  theme_bw()+
  xlab("")+
  coord_flip()+
  ylab("ASV Richness")+
  scale_fill_manual(values = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000'))
