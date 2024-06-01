##Combined leopard frog data
library(ggplot2)
library(vegan)

#read in meta data
meta<-read.delim('PanamaLeopardFrogs/MappingFile_SERDP_Panama.txt', header=T)

#read in OTU table with all leopard frogs
asv_table<-read.delim("PanamaLeopardFrogs/combined_leopard_frog_table.txt", row.names=1, header=T)

#read in OTU table with Panama froggies
asv_table2<-read.delim('PanamaLeopardFrogs/asv_table.txt', header=T, row.names=1)

#read in meta data for Panama froggies
meta2<-read.delim('PanamaLeopardFrogs/filtered_map_Piedra2.txt', header=T)

#filter out non-Rana frogs
meta3<-meta2[which(meta2$Species == "Lithobates warszewitschii"),]
asv_table3<-asv_table2[,which(names(asv_table2) %in% meta3$SampleID)]

#add OTU IDs as new column
asv_table$OTU.ID<-row.names(asv_table)
asv_table3$OTU.ID<-row.names(asv_table3)

#merge all leopard frog data together
MyMerge <- function(x, y){
  df <- merge(x, y, by= "OTU.ID", all.x= TRUE, all.y= TRUE)
  return(df)
}
all_lep_frog <- Reduce(MyMerge, list(asv_table, asv_table3))
all_lep_frog[is.na(all_lep_frog)] <- 0

#write to text file for safe keeping
write.table(all_lep_frog, 'PanamaLeopardFrogs/all_leopard_frog_otu_table.txt', sep='\t', quote=F, row.names=F)

#read in simplified map
meta<-read.delim('PanamaLeopardFrogs/simplified_map_all_samples.txt', header=T)

#calculate beta diversity
ko_pcoa<-capscale(t(all_lep_frog[,-1])  ~ 1, distance='jaccard')

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
#19.2
100*round(ko_pcoa$CA$eig[2]/sum(ko_pcoa$CA$eig), 3)
#6.5

#plot PCoA
ggplot(ko.coords, aes(MDS1, MDS2, colour=Species))+
  geom_point(size=2.9)+
  scale_color_manual(values = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000'))+
  theme_bw()+
  guides(alpha = "none")+
  xlab("PC1- 19.2%")+
  ylab("PC2- 6.5%")

#calculate PERMANOVA with adonis2
row.names(all_lep_frog)<-all_lep_frog$OTU.ID
s16.dis<-as.data.frame(as.matrix(vegdist(t(all_lep_frog[,-1]), method='bray')))
s16.dis2<-s16.dis
s16.dis2$SampleID<-row.names(s16.dis2)

###Calculate pairwise comparisons
#load and reshape data
meta<-read.delim('PanamaLeopardFrogs/simplified_map_all_samples.txt', header=T)
all_lep_frog<-read.delim('PanamaLeopardFrogs/all_leopard_frog_otu_table.txt', header=T, row.names = 1)
s16.dis<-as.data.frame(as.matrix(vegdist(t(all_lep_frog), method = 'jaccard')))

library(reshape2)
s16.dis3<-s16.dis
s16.dis3$variable<-row.names(s16.dis3)
s16.dis_m<-melt(s16.dis3)
names(s16.dis_m)<-c('Sample1', 'Sample2', 'distance')
meta_sub<-meta[,c(1,6)]

#add metadata
s16.dis_m<-merge(s16.dis_m, meta_sub, by.x='Sample1', 'SampleID')
s16.dis_m<-merge(s16.dis_m, meta_sub, by.x='Sample2', 'SampleID')
names(s16.dis_m)<-c('Sample1', 'Sample2', 'Distance', 'Species_2', 'Species_1')

#add new column for comparisons
s16.dis_m$Comparison<-paste(s16.dis_m$Species_2,s16.dis_m$Species_1, sep = ':')

#remove comparisons without Pan lep frog
library(stringr)
lep_frog <- s16.dis_m %>% filter(str_detect(Comparison, 'Bug'))
write.table(lep_frog, 'PanamaLeopardFrogs/pairwise_distance.txt', sep='\t', quote=F, row.names = F)
lep_frog <-read.delim('PanamaLeopardFrogs/pairwise_distance.txt', header=T)

#plot it
ggplot(lep_frog, aes(Comparison, Distance))+
  geom_boxplot()+
  theme_bw()+
  xlab("")+
  ylab('Unweighted UniFrac Similarity')+
  coord_flip()

#bartlett.test(lep_frog$Distance, lep_frog$Comparison)
pairwise.t.test(lep_frog$Distance, lep_frog$Comparison)


#select only columns of interest for the metadata
meta2<-meta[,c(1,5,6,7)]

#add metadata 
s16.dis2<-merge(s16.dis2, meta2, by='SampleID')

#remove sample id column
s16.dis3<-s16.dis2[,2:385]

#calculate adonis
adonis2(s16.dis3[,-c(382,383,384)] ~ Species + Season + Site, data=s16.dis3, permutations = 10000)

#          Df SumOfSqs      R2      F    Pr(>F)    
#Species    6   1.5948 0.37776 42.837 9.999e-05 ***
#Site       5   1.5701 0.37189 45.946 0.009901 **
#Season     3   0.2442 0.05785 13.120 9.999e-05 ***
#Residual 384   2.3827 0.56439                     
#Total    393   4.2218 1.00000 

###########################################################################
#############################Alpha diversity###############################
###########################################################################

#load in asv table
asv.tbl<-read.delim('PanamaLeopardFrogs/all_leopard_frog_otu_table.txt', row.names = 1, header=T)

#read in simplified map
meta<-read.delim('PanamaLeopardFrogs/simplified_map_all_samples.txt', header=T)

#CALCULATE RICHNESS & add metadata & statistics
larv.alph<-as.data.frame(specnumber(t(asv.tbl)))
larv.alph$SampleID<-row.names(larv.alph)
larv.alph<-merge(larv.alph, meta, by='SampleID')
larv.alph$Richness<-as.numeric(larv.alph$`specnumber(t(asv.tbl))`)
pairwise.t.test(larv.alph$Richness, larv.alph$Species, p.adjust.method = 'hochberg')


#plot richness
ggplot(larv.alph, aes(Species, Richness, fill=Species))+
  geom_jitter()+
  geom_boxplot()+
  theme_bw()+
  xlab("")+
  coord_flip()+
  ylab("sOTU Richness")+
  scale_fill_manual(values = c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000'))

#math it
library(broom)
pairwise_result <-pairwise.t.test(larv.alph$`specnumber(t(asv.tbl))`, larv.alph$Species, p.adjust.method = 'BH')

# Convert results to data frame and remove non significant comparisons
pairwise_df1 <- tidy(pairwise_result)
pairwise_df1<-pairwise_df1[-which(pairwise_df1$p.value>0.05),]


###########################################################################
#############################anti-Bd function##############################
###########################################################################
#read in results from vsearch clustering against AmphiBac database
inhibitory<-read.delim("PanamaLeopardFrogs/leopard_frog_inhib_otus.txt", header=F)
inhibitory2<-read.delim("PanamaLeopardFrogs/leopard_frog_out.txt", header=F)
#ASVs are V1 in df

#merge inhibitory datasets
inhibitory_all<-merge(inhibitory, inhibitory2, by='V1', all.x = T, all.y=T)

#load in meta data
meta<-read.delim('PanamaLeopardFrogs/simplified_map_all_samples.txt', header=T)

#load in OTU table
asv.tbl<-read.delim('PanamaLeopardFrogs/all_leopard_frog_otu_table.txt', row.names = 1, header=T)

#subset ASV table to only include ASVs that amtch database
inhib_tb<-asv.tbl[row.names(asv.tbl) %in% inhibitory$V1,]

#calculate percentage of ASVs that are inhibitory
100*(302/15595)
#1.94%

#calculate colSums for total/inhibitory communities
total_sum<-as.data.frame(colSums(asv.tbl))
inhib_sum<-as.data.frame(colSums(inhib_tb))

#bind data together
inhib_tb2<-cbind(total_sum, inhib_sum)

#calculate percent inhibitory
inhib_tb2$per_inhib<-inhib_tb2$`colSums(inhib_tb)`/inhib_tb2$`colSums(asv.tbl)`

#add column for SampleID and merge metadata
inhib_tb2$SampleID<-row.names(inhib_tb2)
inhib_tb2<-merge(inhib_tb2, meta, by='SampleID')

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

#math it
library(broom)
pairwise_result <-pairwise.t.test(inhib_tb2$per_inhib, inhib_tb2$Species, p.adjust.method = 'BH')

# Convert results to dataframe and remove non significant comparisons
pairwise_df <- tidy(pairwise_result)
pairwise_df<-pairwise_df[-which(pairwise_df$p.value>0.05),]
View(pairwise_df)

############################################################################
############################################################################
##compare Bd inhibition to mucosal function for Panama froggies
#read mucosal function data
meta_new<-read.delim('PanamaLeopardFrogs/piedra_map3.txt', header=T)

#read inhibitory data
inhibitory2<-read.delim("PanamaLeopardFrogs/leopard_frog_out.txt", header=F)

#read in OTU table with Panama froggies
asv_table2<-read.delim('PanamaLeopardFrogs/asv_table.txt', header=T, row.names=1)

#remove OTUs that aren't inhibitory
inhib_tb<-asv_table2[row.names(asv_table2) %in% inhibitory2$V1,]

#calculate colSums for total/inhibitory communities
total_sum<-as.data.frame(colSums(asv_table2))
inhib_sum<-as.data.frame(colSums(inhib_tb))

#bind data together
inhib_tb2<-cbind(total_sum, inhib_sum)

#calculate percent inhibitory
inhib_tb2$per_inhib<-inhib_tb2$`colSums(inhib_tb)`/inhib_tb2$`colSums(asv_table2)`

#add column for SampleID and merge metadata
inhib_tb2$SampleID<-row.names(inhib_tb2)
inhib_tb2<-merge(inhib_tb2, meta_new, by='SampleID')

#remove NA's for mucosal function
inhib_tb2<-inhib_tb2[complete.cases(inhib_tb2$Mucosome.Bd.inhibition.per.cm2.skin),]

#calculate correlations bewteen Bd inhibition of mucosome and Bd inhibition percent
library(lme4)

muc_model<-lmer(inhib_tb2$Mucosome.Bd.inhibition.per.cm2.skin ~ inhib_tb2$per_inhib | inhib_tb2$Species)
summary(muc_model)
anova(muc_model)

ggplot(inhib_tb2, aes(Mucosome.Bd.inhibition.per.cm2.skin, per_inhib))+
         theme_bw()+
         geom_point()+
  geom_smooth(method = 'lm')+
  ylab("Percent Inhibitory Towards Bd")+
  xlab("Mucosal function (Bd inhibition cm^2)")


###########################################################################
#################Percent occupancy#########################################
###########################################################################
library(ggplot2)
library(plyr)
library(reshape2)
library(vegan)

#load in OTU table
asv.tbl<-read.delim('PanamaLeopardFrogs/all_leopard_frog_otu_table.txt', row.names = 1, header=T)

#read in simplified map
meta<-read.delim('PanamaLeopardFrogs/simplified_map_all_samples.txt', header=T)

#calculate log OTU abundance
otu_means<-rowMeans(asv.tbl)
otu_means<-as.data.frame(otu_means)

#take log of means
otu_means_log<-log(rowMeans(asv.tbl))
otu_means_log<-as.data.frame(otu_means_log)

#calculate the presence of OTUs
otu_table_binary<-asv.tbl
otu_table_binary[otu_table_binary>0]<-1
otu_samp_count<-rowSums(otu_table_binary)
otu_samp_count<-as.data.frame(otu_samp_count)

#calculate percent occupancy
dim(asv.tbl)
#15595 OTUs, 381 samples

#make a table of OTU counts and mean abundance
per_occu<-cbind(otu_samp_count, otu_means, otu_means_log)
per_occu$per_occ<-((100*(per_occu$otu_samp_count/164)))

#remove -Inf
per_occu <- subset(per_occu, !is.infinite(otu_means_log))
per_occu$occup<-per_occu$otu_samp_count/381

#plot all
ggplot(per_occu, aes(occup, otu_means_log))+
  geom_point()+
  theme_bw()+
  ylab('Log sOTU Abundance')+
  xlab('Percent Occupancy')


#get OTUs of interest (>50% of sample and log abundance >2)
count55<-which(per_occu$occup>0.5)
count55_p<-per_occu[count55,]

count3<-which(per_occu$otu_means_log>2)
count3_p<-per_occu[count3,]

ggplot(per_occu, aes(occup, otu_means_log))+
  geom_point(aes(colour='grey', size=2))+
  geom_point(count55_p, aes(occup, otu_means_log, size=2), colour='black')+
  geom_point(count3_p, aes(occup, otu_means_log, size=2), colour='black')+
  theme_bw()+
  ylab('Log sOTU Abundance')+
  xlab('Percent Occupancy')

#make a data frame of 'core' OTUs
core_otus<-rbind(count3_p, count55_p)
core_otus$OTU<-row.names(core_otus)
dim(core_otus)
#[1] 66  5

#add taxonomy to see who they are
tax<-read.delim('PanamaLeopardFrogs/all_leopard_frog_taxonomy.tsv')

#bind taxonomy
core_otus<-merge(core_otus, tax, by.x = 'OTU', by.y='Feature.ID', all.x=T, all.y=F)
core_otus<-na.omit(core_otus[!is.na(core_otus$Taxon), ])

#add inhibitory data
inhibitory<-read.delim("PanamaLeopardFrogs/leopard_frog_inhib_otus.txt", header=F)
inhibitory2<-read.delim("PanamaLeopardFrogs/leopard_frog_out.txt", header=F)
#ASVs are V1 in df

#merge inhibitory datasets
inhibitory_all<-merge(inhibitory, inhibitory2, by='V1', all.x = T, all.y=T)

#add info to core OTU table
core_otus<-merge(core_otus, inhibitory_all, all.x = T, all.y=F, by.x='OTU', by.y='V1')

#add abundance/counts for each sample
library(dplyr)
library(plyr)
library(reshape2)

asv2<-asv.tbl
asv2$OTU<-row.names(asv2)
asv_m<-melt(asv2)
asv_m<-asv_m[asv_m$OTU %in% core_otus$OTU,]
asv_m<-merge(asv_m, meta, by.x='variable', by.y='SampleID')
asv_m$value<-as.numeric(asv_m$value)

asv_sum<-ddply(asv_m, c('OTU', 'Species'), summarize, mean=mean(value), n=100*length(which((value>=1)))/length(value))
asv_sum2<-dcast(asv_sum, formula = OTU ~ Species)

#add to core OTU data
core_otus2<-merge(core_otus, asv_sum2, by.y='OTU', by.x='OTU')

#write to file
write.table(core_otus2, 'PanamaLeopardFrogs/core_lep_frog_otus.txt', quote=F, sep='\t', row.names=F)

#######################################################
#######################Taxa plot#######################
#######################################################
taxa_data<-read.delim('PanamaLeopardFrogs/new axa_plot_data.txt')
library(ggplot2)
library(reshape2)

#reshape data for plotting
taxa_m<-melt(taxa_data)
names(taxa_m)<-c('Species', 'Order', 'value')

#read in the best pallette
pal<-c("#771155", "#CC99BB", "#114477", "#4477AA", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788","#41AB5D", "#252525", "#525252", "#737373", "#969696")

#plot data
ggplot(taxa_m, aes(Species, value, fill=Order))+
  geom_bar(stat='identity')+
  scale_fill_manual(values=pal)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  ylab("Relative Abundance")+
  xlab("")
