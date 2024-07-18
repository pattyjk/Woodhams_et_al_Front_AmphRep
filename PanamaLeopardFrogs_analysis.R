#load packages
library(vegan)
library(ggplot2)
library(ggpubr)
library(plyr)
library(dplyr)
library(reshape2)
library(RColorBrewer)
###########################################################################
#############################Beta diversity################################
###########################################################################

#load in OTU table
asv.tbl<-read.delim('PanamaLeopardFrogs/asv_table.txt', row.names = 1, header=T)

#load in meta data
meta<-read.delim('PanamaLeopardFrogs/filtered_map_Piedra2.txt', header=T)

#subset asv table to remove non-Alto de Piedra samples (that's all that's in the filtered map)
asv.tbl<-asv.tbl[,names(asv.tbl) %in% meta$SampleID]

#get column sums for rarefaction
min(colSums(asv.tbl))
#1322

#rarefy table
set.seed(515)
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
ggplot(ko.coords, aes(MDS1, MDS2, colour=Species, size=Log_Bd))+
  geom_point()+
  scale_color_manual(values = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000'))+
  theme_bw()+
  guides(alpha = "none")+
  xlab("PC1- 32.2%")+
  ylab("PC2- 12.6%")

#plot axis1 against log Bd load
ggplot(ko.coords, aes(MDS1,  Log_Bd))+
  geom_point()+
  theme_bw()+
 stat_cor(method = "spearman", cor.coef.name="rho")+
  geom_smooth(method = 'lm')+
  ylab("Log Bd load")+
  xlab("Axis 1 values")
#Significant negative relationship, rho=-0.55, p=6.7e-7

ggplot(ko.coords, aes(MDS2, Log_Bd))+
  geom_point()+
  theme_bw()+
  stat_cor(method = "spearman", cor.coef.name="rho")+
  geom_smooth(method = 'lm')+
  ylab("Bd load")+
  xlab("Axis 2 values")
#significant rho=-0.26, p=0.015
  
#calculate non-parametric permanova
#calculate distance
s16.dis<-as.data.frame(as.matrix(vegdist(asv.rare, method='jaccard')))
s16.dis2<-s16.dis
s16.dis2$SampleID<-row.names(s16.dis2)

#select only columns of interest for the metadata
meta2<-meta[,c(1,9,17)]

#add metadata 
s16.dis2<-merge(s16.dis2, meta2, by='SampleID')
  
#remove sample id column
s16.dis3<-s16.dis2[,2:92]
  
#calculate adonis
adonis2(s16.dis3[,-c(90,91)] ~ Species + Bd.GE, data=s16.dis3, permutations = 10000)

#          Df SumOfSqs      R2      F    Pr(>F)    
#Species   9  0.19297 0.38741 5.5511 9.999e-05 ***
#Bd.GE     1  0.0211 0.1419 1.5353    0.007437 
#Residual 79  0.30515 0.61259                     
#Total    88  0.49812 1.00000         


###########################################################################
#############################Alpha diversity###############################
###########################################################################

#load in asv table
asv.tbl<-read.delim('PanamaLeopardFrogs/asv_table.txt', row.names = 1, header=T)

#subset asv table to remove non-Alo de Piedra samples (that's all that's in the filtered map)
asv.tbl<-asv.tbl[,names(asv.tbl) %in% meta$SampleID]

#load in meta data
meta<-read.delim('PanamaLeopardFrogs/filtered_map_Piedra2.txt', header=T)

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
  ylab("sOTU Richness")+
  scale_fill_manual(values = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000'))

bartlett.test(larv.alph$`specnumber(t(asv.tbl))`, larv.alph$Species)
#Bartlett's K-squared = 17.196, df = 9, p-value = 0.04573

pairwise.t.test(larv.alph$`specnumber(t(asv.tbl))`, larv.alph$Species, p.adjust.method = 'hochberg')
#no sig differences

###########################################################################
#############################anti-Bd function##############################
###########################################################################
#read in results from vsearch clustering against AmphiBac database
inhibitory<-read.delim("PanamaLeopardFrogs/leopard_frog_inhib_otus.txt", header=F)
#ASVs are V1 in df

#load in meta data
meta<-read.delim('PanamaLeopardFrogs/piedra_map3.txt', header=T)

#load in asv table
asv.tbl<-read.delim('PanamaLeopardFrogs/asv_table.txt', row.names = 1, header=T)

#subset asv table to remove non-Alo de Piedra samples (that's all that's in the filtered map)
asv.tbl<-asv.tbl[,names(asv.tbl) %in% meta$SampleID]

#subset ASV table to only include ASVs that amtch database
inhib_tb<-asv.tbl[row.names(asv.tbl) %in% inhibitory$V1,]

#CALCULATE RICHNESS & add metadata & statistics
library(vegan)
larv.alph2<-as.data.frame(specnumber(t(inhib_tb)))
larv.alph2$SampleID<-row.names(larv.alph2)
larv.alph2<-merge(larv.alph2, meta, by='SampleID')
larv.alph2$Richness<-as.numeric(larv.alph2$`specnumber(t(inhib_tb))`)

#plot it
library(ggplot2)
ggplot(larv.alph2, aes(Mucosome.Bd.inhibition.per.cm2.skin, Richness))+
  geom_point()+
  theme_bw()+
  ylab("Antifungal Richness")+
  ylab("Mucosome function (per cm2 skin)")

#extract onlt Pan lep frog
pan.lep<-larv.alph2[which(larv.alph2$Species == "Ngäbe-Buglé leopard frog"),]
ggplot(pan.lep, aes(Mucosome.Bd.inhibition.per.cm2.skin, Richness))+
  geom_point()+
  theme_bw()+
  ylab("Antifungal Richness")+
  ylab("Mucosome function (per cm2 skin)")

cor.test(pan.lep$Richness, pan.lep$Mucosome.Bd.inhibition.per.cm2.skin)
cor.test(larv.alph2$Richness, larv.alph2$Mucosome.Bd.inhibition.per.cm2.skin)

#calculate colSums for total/inhibitory communities
total_sum<-as.data.frame(colSums(asv.tbl))
inhib_sum<-as.data.frame(colSums(inhib_tb))

#bind data together
inhib_tb2<-cbind(total_sum, inhib_sum)

#calculate percent inhibitory
inhib_tb2$per_inhib<-inhib_tb2$`colSums(inhib_tb)`/inhib_tb2$`colSums(asv.tbl)`

#change to out of 100%
inhib_tb2$per_inhib<-inhib_tb2$per_inhib*100

#add column for SampleID and merge metadata
inhib_tb2$SampleID<-row.names(inhib_tb2)
inhib_tb2<-merge(inhib_tb2, meta, by='SampleID')

#merge inhibi with richness
inhib_tb2<-merge(inhib_tb2, larv.alph2, by='SampleID')



#plot per species
library(ggplot2)
ggplot(inhib_tb2, aes(Species, per_inhib, fill=Species))+
  geom_jitter()+
  geom_boxplot()+
  theme_bw()+
  xlab("")+
  coord_flip()+
  ylab("Percent Inhibitory towards Bd")+
  scale_fill_manual(values = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000'))

ggplot(inhib_tb2, aes(Mucosome.Bd.inhibition.per.cm2.skin, per_inhib))+
  geom_point()+
  theme_bw()+
  ylab("Percent inhibitory towards Bd")+
  ylab("Mucosome function (per cm2 skin)")

cor.test(inhib_tb2$per_inhib, inhib_tb2$Mucosome.Bd.inhibition.per.cm2.skin)

#calculate stats
pairwise.t.test(inhib_tb2$per_inhib, inhib_tb2$Species, p.adjust.method = 'hochberg')
#most comparisions are ns, except HyCo vs SmSi, Ngäbe-Buglé leopard frog vs. SmSi, LiWa cs SmSi


#################Look at taxa that respond to Bd load
#load in OTU table
asv.tbl<-read.delim('PanamaLeopardFrogs/asv_table.txt', row.names = 1, header=T)

#load in meta data
meta<-read.delim('PanamaLeopardFrogs/filtered_map_Piedra2.txt', header=T)

#subset asv table to remove non-Alto de Piedra samples (that's all that's in the filtered map)
asv.tbl<-asv.tbl[,names(asv.tbl) %in% meta$SampleID]

#get column sums for rarefaction
min(colSums(asv.tbl))
#1322

#rarefy table
set.seed(515)
asv.rare<-rrarefy(t(asv.tbl), sample=1322)

#reshape OTU table
asv_m<-melt(asv.rare)

#add meta data
asv_m<-merge(asv_m, meta, by.x='Var1', by.y='SampleID', all.y = F, all.x=T)

#calculate correlations
library(dplyr)

calculate_spearman <- function(asv_m) {
  cor_test <- cor.test(asv_m$Bd.GE, asv_m$value, method = "spearman")
  return(data.frame(correlation = cor_test$estimate, p_value = cor_test$p.value))
}

results <- asv_m %>%
  group_by(Var2) %>%
  do(calculate_spearman(.))

#correct p-value with Hochberg
results$p_corr<-p.adjust(results$p_value, method = 'hochberg')

#how many
length(which(results$p_value<0.05))
#148
length(which(results$p_corr<0.05))
#14

#create table that only has significant ones
sig_otus<-results[which(results$p_corr<0.05),]

#add taxonomy
tax<-read.delim("PanamaLeopardFrogs/all_leopard_frog_taxonomy.tsv")
sig_otus<-merge(sig_otus, tax, by.x='Var2', by.y='Feature.ID')

#inhibitory towards Bd?
inhib_otus<-read.delim('PanamaLeopardFrogs/leopard_frog_inhib_otus.txt', header=F)
sig_otus<-merge(sig_otus, inhib_otus, by.x='Var2', by.y='V1', all.x=T, all.y=F)

#write to table
write.table(sig_otus, 'PanamaLeopardFrogs/otus_bd_load.txt', sep='\t', quote=F, row.names = F)

