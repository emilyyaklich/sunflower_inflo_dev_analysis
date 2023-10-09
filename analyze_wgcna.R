# Name: analyze wgcna
# Author: EY
# Date: 06/29/2023
# Version:4.2.1
# Description: will analyze wgcna network for the cultivated expression analysis....ANOVA?

#install.packages("BiocManager")

library(WGCNA)
library(mltools)
library(data.table)
library(FDRestimation)
library(car)
library(emmeans)
library(gtools)
library(dplyr)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(ggpubr)
library(goseq)
library(GO.db)
library(tidyr)

setwd('/home/ely67071/sunflower_stress_analysis/')


# read in the cultivated metadata for wgcna
metadata<-read.csv(file='wgcna/combined_metadata_wgcna.csv')

input_mat_filt<-readRDS("wgcna/input_matrix.RData")

table(dimnames(input_mat_filt)[[1]]==metadata$Plant.)

# read in the network
netwk<-readRDS("wgcna/netwk_bicor_blocksize.RData")

# columns and rows of the matrix
nGenes<-ncol(input_mat_filt)
nsamples<-nrow(input_mat_filt)
# turn the labels into colors
mergedColors=labels2colors(netwk$colors)
netwk$colors


# create a table of module sizes
module_sizes <- as.data.frame(table(mergedColors))
module_sizes$mergedColors <-paste("ME", module_sizes$mergedColors,sep="")


MEs0 <- moduleEigengenes(input_mat_filt, mergedColors)$eigengenes
MEs <- orderMEs(MEs0)
help("orderMEs")
treatments<-metadata[,c(11)]


treatments<-as.factor(treatments)
treatments_onehot<-one_hot(as.data.table(treatments))

all_treatments<-as.integer(grepl(pattern="Control",x=treatments))
all_treatments<-1-all_treatments

treatments_onehot<-cbind(treatments_onehot,all_treatments)

colnames(treatments_onehot)
colnames(treatments_onehot)<-c("Combo", "Control","HighSalt","LowNut", "AllStress")


treatments_onehot<-treatments_onehot[,c("Control","AllStress","HighSalt", "Combo","LowNut")]

# remove control and all stress
#treatments_onehot<-treatments_onehot[,-c(1:2)]
#treatments_onehot_combo_split<-treatments_onehot
# use this loop to 
#for(row in 1:nrow(treatments_onehot_combo_split)){
#  if (treatments_onehot_combo_split[row]$Combo==1){
#    treatments_onehot_combo_split[row]$HighSalt <-1
#  }
#  if (treatments_onehot_combo_split[row]$Combo==1){
#    treatments_onehot_combo_split[row]$LowNut <-1
#  }
#}
#treatments_onehot_combo_split<-treatments_onehot_combo_split[,-c(2)]


moduleTraitCor<-cor(MEs, treatments_onehot, use="p")
moduleTraitpval<-corPvalueStudent(moduleTraitCor, nsamples)

#png('wgcna/trait-correlations.png', width=2000, height =2200, res=300)

#sizeGrWindow(6,10)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitpval, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(treatments_onehot),
              yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1))
               #main = paste("Module-trait relationships"))
#dev.off()



# linear regression 
# process the metadata
metadata$Treatment<-relevel(as.factor(metadata$Treatment),ref="Control")

levels(metadata$Treatment)
metadata$Accession<-as.factor(make.names(metadata$Accession))
levels(metadata$Accession)

metadata$Type<-as.factor((metadata$Type))
levels(metadata$Type)

metadata$SampleDay<-factor((metadata$SampleDay))
levels(metadata$SampleDay)

metadata$Bench<-factor(metadata$Bench)
levels(metadata$Bench)

metadata$Reproductive<-factor(metadata$Reproductive)
levels(metadata$Reproductive)

# samples per accession/treatment
aggregate(metadata$Plant., by=list(metadata$Accession, metadata$Treatment),length)

grep("Type", colnames(metadata))
grep("Accession", colnames(metadata))
grep("Salt", colnames(metadata))
grep("Osmocote", colnames(metadata))
grep("Treatment", colnames(metadata))
grep("Reproductive", colnames(metadata))
grep("SampleDay", colnames(metadata))
grep("Bench", colnames(metadata))


metadata$Osmocote<-relevel(as.factor(metadata$Osmocote),ref="High")
metadata$Salt<-relevel(as.factor(metadata$Salt),ref="NoSalt")



all <-metadata[,c(5,77,4,10,9,11,63,64,6)]

MEs$"Plant." <-rownames(MEs)
MEs$Plant.

AllData<-merge(metadata[,c(5,77,4,10,9,11,63,64,6)], MEs, by="Plant.")


levels(AllData$Osmocote)
levels(AllData$Salt)

# think about recoding this
LR_Mod <- function(Yvar,dataset) {
  mod <- lm(Yvar ~Osmocote + Salt + Osmocote:Salt + Type + Type:Accession + 
              Osmocote:Type + Salt:Type + Osmocote:Salt:Type + Bench,
            data=dataset, contrast=list(Accession=contr.sum, Bench=contr.sum))
  return(mod)
}


LR_mod_results <- lapply(AllData[,c(10:59)], function(x) {LR_Mod(x,AllData)})

str(LR_mod_results)



LR_mod_ANOVA <- lapply(LR_mod_results, function(x) {as.data.frame(Anova(x,test="F",type=2))})

# save f and p vals
anova_columns <- lapply(LR_mod_ANOVA, function(x) {x[c(1:3,5,6),c(3,4)]})



anova_wide<-lapply(anova_columns,function(x){as.data.frame(t(x))})
anova_widewlabels<-lapply(names(anova_wide),function(x) {anova_wide[[x]]$Module <- x;return(anova_wide[[x]])} ) 

#f vals
anova_fvals <- lapply(anova_widewlabels, function(x) {x[1,]})

# combine into 1 df
all_fvals <- do.call("rbind", anova_fvals)
colnames(all_fvals) <- c("nut_f","salt_f","type_f","nutxsalt_f", "typexaccession_f","Module")

anova_pvals <- lapply(anova_widewlabels, function(x) {x[2,]})
all_pvals <- do.call("rbind", anova_pvals)
colnames(all_pvals) <- c("nut_p","salt_p","type_p","nutxsalt_p", "typexaccession_p","Module")

anova_results <- merge(all_fvals, all_pvals, by="Module")
write.csv(anova_results, file="wgcna/wgcna_anova_results.csv", row.names=FALSE)


# LS means




# test osmocote and salt differently 
#mod_means_osmocote<-lapply(LR_mod_results, function(x) {emmeans(x, ~Osmocote, type="response")})
#mod_means_salt<-lapply(LR_mod_results, function(x) {emmeans(x, ~Salt, type="response")})

#mod_means_df_salt<-lapply(mod_means_salt, function(x) {as.data.frame(x)[,c(1:4)]})
#mod_means_w_labels_salt<-lapply(names(mod_means_df_salt), function(x) {mod_means_df_salt[[x]]$Module <- x;return(mod_means_df_salt[[x]])})
#all_mod_means_salt<-do.call("rbind", mod_means_w_labels_salt)

#all_mod_means_salt$Treatment <- ifelse(all_mod_means_salt$Salt == "NoSalt" ,"Control_Salt", "Salt")


#mod_means_df_osmocote<-lapply(mod_means_osmocote, function(x) {as.data.frame(x)[,c(1:4)]})
#mod_means_w_labels_osmocote<-lapply(names(mod_means_df_osmocote), function(x) {mod_means_df_osmocote[[x]]$Module <- x;return(mod_means_df_osmocote[[x]])})
#all_mod_means_osmocote<-do.call("rbind", mod_means_w_labels_osmocote)

#all_mod_means_osmocote$Treatment <- ifelse(all_mod_means_osmocote$Osmocote == "High" ,"Control_Nut", "LowNut")

#test<-merge(all_mod_means_osmocote, all_mod_means_salt, by="Module")
#all_mod_means_wide_osmocote <-reshape(all_mod_means_osmocote[,c(2:3,5:6)], idvar="Module", timevar = "Treatment", direction="wide")
#all_mod_means_wide_salt <-reshape(all_mod_means_salt[,c(2:3,5:6)], idvar="Module", timevar = "Treatment", direction="wide")
#all_mod_means_wide_salt_osmocote<-merge(all_mod_means_wide_osmocote, all_mod_means_wide_salt, by="Module")


# get a note for this that nested model....



# below is if we want to use combo in the model!!

# for all modules
mod_means<-lapply(LR_mod_results, function(x) {emmeans(x, ~Osmocote*Salt, type="response")})
mod_means_test<-lapply(LR_mod_results, function(x) {emmeans(x, ~Osmocote*Salt, type="response")})

#identical(mod_means,mod_means_test)
mod_means_df<-lapply(mod_means, function(x) {as.data.frame(x)[,c(1:4)]})

mod_means_w_labels<-lapply(names(mod_means_df), function(x) {mod_means_df[[x]]$Module <- x;return(mod_means_df[[x]])})
all_mod_means<-do.call("rbind", mod_means_w_labels)

all_mod_means$Treatment <- ifelse(all_mod_means$Osmocote == "High" & all_mod_means$Salt == "NoSalt", "Control", ifelse(all_mod_means$Osmocote=="High" & all_mod_means$Salt == "Salt", "Salt", ifelse(all_mod_means$Osmocote =="Low" & all_mod_means$Salt =="NoSalt", "LowNut", "Combo")))


# check
aggregate(all_mod_means$Module, by=list(all_mod_means$Osmocote, all_mod_means$Salt, all_mod_means$Treatment),length)

all_mod_means_wide <-reshape(all_mod_means[,c(3:6)], idvar="Module", timevar = "Treatment", direction="wide")




# combine and save results

all_results<-merge(anova_results, all_mod_means_wide, by="Module")


# adjust p-value using Bonferroni-Holm (also known as seqeuential bonferroni)
all_pvals_adj<-lapply((all_results[7:11]), function(x) {p.adjust(x, method="holm", n=length(x))})
all_pvals_adj<-as.data.frame(do.call(cbind, all_pvals_adj))

# rename columns
colnames(all_pvals_adj) <- c("nut_p_adj","salt_p_adj","type_p_adj","nutxsalt_p_adj", "typexaccession_p_adj")

# combine the dataframes
all_results_adj<- cbind(all_results, all_pvals_adj)

# turn the p-values into significance stars
all_stars<-lapply((all_results_adj[20:24]), function(x) {stars.pval(x)})
all_stars<-as.data.frame(do.call(cbind, all_stars))

# name columns
colnames(all_stars) <- c("nut_stars","salt_stars","type_stars","nutxsalt_stars", "typexaccession_stars")
all_results_stars <- cbind(all_results_adj, all_stars)

# set the stars as factors
all_results_stars$nut_stars<-factor(all_results_stars$nut_stars, levels=c("***", "**",   "*",  " " ))
all_results_stars$salt_stars<-factor(all_results_stars$salt_stars, levels=c("***", "**",   "*",  " " ))
all_results_stars$nutxsalt_stars<-factor(all_results_stars$nutxsalt_stars, levels=c("***", "**",   "*",  " " ))

# remove the NA
all_results_stars[is.na(all_results_stars)] <- " "

# sort the columns based on signficance (star levels)
all_results_stars_sorted<-all_results_stars %>% arrange(nut_stars, salt_stars, nutxsalt_stars)

# drop the type and typexaccession
drops<-c("type_stars","typexaccession_stars")
all_results_stars_sorted<-all_results_stars_sorted[,!names(all_results_stars_sorted) %in% drops]

# subset the columns we wnat to plot
plotting_subset<-all_results_stars_sorted[,c('Module',"emmean.Control","SE.Control","emmean.LowNut","SE.LowNut","emmean.Salt", "SE.Salt", "emmean.Combo", "SE.Combo","emmean.Salt", "nut_stars", "salt_stars", "nutxsalt_stars")]

# extract module names from the columns
module_names<-(plotting_subset[,1])




# loop through and create variables to plot each module
plot_list_nvs<-list()
i<-0
for(module in module_names){
  i <- i+1
  module_subset<-subset(plotting_subset, Module %in% module)
  
  plotting_df_wrt_nutrient<-data.frame(Nutrient=c("Control", "Stress"), salt_control=c(module_subset$emmean.Control,module_subset$emmean.LowNut),salt_stress=c(module_subset$emmean.Salt,module_subset$emmean.Combo))
  data_long<-melt(plotting_df_wrt_nutrient, id="Nutrient")
  colnames(data_long) <- c("Nutrient", "Stress_type", "emmean_eigenvector")
  
  module_size=subset(module_sizes, mergedColors ==module)
  plot_title<-paste(module_subset$Module, module_size$Freq, sep=": ")
  
  
  module_size=subset(module_sizes, mergedColors ==module)
  plot_title<-paste(module_subset$Module, module_size$Freq, sep=": ")
  
  nut<-paste("Nut: ", as.character(module_subset$nut_stars))
  salt<-paste("Salt: ", as.character(module_subset$salt_stars))
  nutxsalt<-paste("NxS: ", as.character(module_subset$nutxsalt_stars))
  
  
  
  
  
  plot_list_nvs[[i]]<-ggplot(data_long, aes(x=Nutrient, y=emmean_eigenvector, color=Stress_type, group=Stress_type)) +geom_line(size=1)+ scale_color_manual(values=c('blue','darkorange3')) +
    ggtitle(plot_title) + annotate('text',x="Stress",y=0.2,label=nut, size =4) + annotate('text',x="Stress",y=0.16,label=salt, size=4) + annotate('text',x="Stress",y=0.11,label=nutxsalt, size=4) + ylim(-0.2,0.2)+
    theme_bw() + theme(plot.title = element_text(size=20),axis.line=element_line(colour = "black"),
                       panel.grid.major=element_blank(),
                       panel.grid.minor = element_blank(),
                       panel.border=element_blank(),
                       panel.background = element_blank())
  
  path_to_plots1<-paste('wgcna/interaction_plots/individaul_plots/',module_subset$Module,sep="")
  path_to_plots2<-paste(path_to_plots1,".png",sep="")
  png(path_to_plots2, width=100, height =100, res=300)
  ggplot(data_long, aes(x=Nutrient, y=emmean_eigenvector, color=Stress_type, group=Stress_type)) +geom_line(size=1)+ scale_color_manual(values=c('blue','darkorange3')) +
    ggtitle(plot_title) + annotate('text',x="Stress",y=0.2,label=nut, size =4) + annotate('text',x="Stress",y=0.16,label=salt, size=4) + annotate('text',x="Stress",y=0.11,label=nutxsalt, size=4) + ylim(-0.2,0.2)+
    theme_bw() + theme(plot.title = element_text(size=20),axis.line=element_line(colour = "black"),
                       panel.grid.major=element_blank(),
                       panel.grid.minor = element_blank(),
                       panel.border=element_blank(),
                       panel.background = element_blank())
  dev.off()
}

# print individual plots
j<-0
for (module in module_names) {
  j <- j+1
  path_to_plots<-paste('wgcna/interaction_plots/individual_plots/',module,".png", sep="")
  png(path_to_plots, width=1000, height =1000, res=300)
  print(plot_list_nvs[[j]])
  dev.off()
}

# plot them in subsets (50 plots seems way too large)
subset_1_nvs<-plot_list_nvs[1:12]
png('wgcna/interaction_plots/subset1_interaxn.png', width=2000, height =2200, res=300)
plot_1_nvs <- do.call(grid.arrange, subset_1_nvs)
dev.off()


subset_2_nvs <- plot_list_nvs[13:24]
png('wgcna/interaction_plots/subset2_interaxn.png', width=2000, height =2200, res=300)
plot_2_nvs <- do.call(grid.arrange, subset_2_nvs)
dev.off()

subset_3_nvs<- plot_list_nvs[25:36]
png('wgcna/interaction_plots/subset3_interaxn.png', width=2000, height =2200, res=300)
plot_3_nvs<-do.call(grid.arrange, subset_3_nvs)
dev.off()

subset_4_nvs<- plot_list_nvs[37:48]
png('wgcna/interaction_plots/subset4_interaxn.png', width=2000, height =2200, res=300)
plot_4_nvs<-do.call(grid.arrange, subset_4_nvs)
dev.off()

subset_5_nvs <- plot_list_nvs[49:50]
png('wgcna/interaction_plots/subset5_interaxn.png', width=2000, height =2200, res=300)
plot_5_nvs<-do.call(grid.arrange, subset_5_nvs)
dev.off()


# plot by result category
print(module_names)



















# GO analysis

# match colors to gene name
color2gene = data.frame(unlist(colnames(input_mat_filt)),unlist(mergedColors))



ha412_genes<-read.table("Ha412GO_terms_interpro_aa_slim.txt", fill=T)
colnames(ha412_genes)<-c("parent", "ontology_term")

length_test <- ha412_genes$parent
length_test2 <- unique(length_test)


ha412_genes$ontology_term<-as.character(ha412_genes$ontology_term)
go_mapping<-separate_rows(ha412_genes, ontology_term, sep=",")
unique_go_mapping<-go_mapping %>% distinct()
distinct(go_mapping)
?distinct()
probes = colnames(input_mat_filt)

test <- ha412_genes %>% dplyr::slice(match(ha412_genes$parent, color2gene$unlist.colnames.input_mat_filt..))
test2<-unique(test)



test3<-inner_join(test2, color2gene, by=c('parent'='unlist.colnames.input_mat_filt..'))

test3<-separate_rows(test3, ontology_term, sep=",")
unique_go_test3<-test3 %>% distinct()
test31 <-distinct(test3)
test_list<-test3$parent
test_list_unique<-unique(test_list)

write.csv(unique_go_test3, "wgcna/go_module_data.csv")



# run GSEA













ha412_genes$parent
module_names
probes2annot = match(ha412_genes$parent, probes)
probes2annot
all_ids = ha412_genes$ontology_term[probes2annot]

netwk$colors

for (module in module_names){
  mod_genes = (moduleColors==module)
  mod_ids = all_ids[mod_genes]
  #print(module)
  print(mod_genes)
}
all_ids

length(mergedColors)
length(all_ids)

GOenr = GOenrichmentAnalysis(mergedColors, all_ids)
