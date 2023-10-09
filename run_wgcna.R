# Name: wgcna
# Author: EY
# Date: 06/20/2023
# Version:4.2.1
# Description: will run wgcna analysis with the cultivated expression data (HA and RHA only)

setwd('/home/ely67071/sunflower_inflo_dev_analysis/')

install.packages("WGCNA")
library(DESeq2)
library(WGCNA)
library(tidyverse)
library(BiocParallel)


samples=c("10D_REP1_ATTACTCG", "20D_REP2_TCCGGAGA" ,"30D_REP2_CGCTCATT", "35D_REP1_GAGATTCC", 
          "HA_10D_2_ACCTTGGC", "HA_10D_3_ATATCTCG", "HA_20D_2_GCGCTCTA", 
          "HA_20D_3_AACAGGTT", "HA_30D_2_GGTGAACC", "HA_30D_3_CAACAATG", "HA_35D_2_TGGTGGCA", "HA_35D_3_AGGCAGAG")

dev_stage<-sub(".*([0-9]{2,2}D).*", "\\1",samples)
dev_stage <- as.factor(dev_stage)

metadata<-data.frame(samples, dev_stage)


# read in deseq results (output from run_DGE_deseq_cult_only.R)
# note this is reading in all of the data and is NOT filtered for genes that are
# differentially expressed 
deseq <- readRDS('deseq_results/combatseq/pairwise/deseq_dataset_results_combatseq.RData')
deseq$samples
#deseq$samples <- paste(deseq$dev_stage, "-", deseq$samples)
#treatment_labels<-deseq$samples


# pre-filter for reads where at least 10 samples have a count of 10 or higher
keep<-rowSums(counts(deseq)>=10)>=10
length(which(keep==1))
deseq_filt<-deseq[keep,]

# transform data into a matrix
vsd <- vst(deseq_filt,blind=TRUE)
vsd_matrix<-assay(vsd)

input_mat=t(vsd_matrix)
input_mat[1:5,1:10]

# check samples
gsg<-goodSamplesGenes(input_mat,verbose=3)
gsg$allOK # this statement gives TRUE, so no further removal is necessary


# see if there are outliers to cut in teh dendrogram
sampleTree<-hclust(dist(input_mat), method="average")
par(cex=0.6)
par(mar=c(0,4,2,0))
png("wgcna/sample_clustering_outliers.png")
plot(sampleTree, main= "Sample clustering to detect outliers", sub="",xlab="",cex.lab=1.5,cex.axis=1.5,cex.main=2)
abline(h=200, col="red")
dev.off()

# remove the outlier samples
clust<-cutreeStatic(sampleTree,cutHeight = 200,minSize=10)
table(clust)
keepSamples<-(clust==1)
input_mat_filt<-input_mat[keepSamples,]
saveRDS(input_mat_filt,file='wgcna/input_matrix.RData')


# remove the outliers from the metadata as well
#remove_outlier<-setdiff(rownames(input_mat),rownames(input_mat_filt))
#metadata<- metadata[!metadata$Plant. %in% remove_outlier,]
#write.csv(metadata, "wgcna/metadata_wgcna.csv")


# get traits of interest from metadata (aka Treatment)
samples<-rownames(input_mat_filt)
trait_rows<-match(samples, metadata$samples)
datTraits<-metadata[trait_rows,-1]
rownames(datTraits)<-metadata[trait_rows,1]
metadata[trait_rows,1]
metadata[trait_rows,-1]
# create treatment labels based after the outliers have been removed
treatment_labels<-paste(metadata$Treatment,"-",metadata$Plant.)

# see if there are global differences between samples 
sampleTree_filt<-hclust(dist(input_mat_filt), method="average")
sample_names<-rownames(input_mat_filt)

# based on stage
traitcolors<-labels2colors(metadata$dev_stage)
png("wgcna/sample_relatedness_clusters.png")
plotDendroAndColors(sampleTree_filt,traitcolors)
# no strong clusters indicating globally different groups of samples 
dev.off()

# based on 
traitcolors2<-labels2colors(metadata$Accession)
png("wgcna/sample_relatedness_clusters_accession.png")
plotDendroAndColors(sampleTree_filt,traitcolors2)
dev.off()

# based on type (heterotic group)
traitcolors3<-labels2colors(metadata$Type)
png("wgcna/sample_relatedness_clusters_type.png")
plotDendroAndColors(sampleTree_filt,traitcolors3)
dev.off()

# ID power
#allowWGCNAThreads()
powers=c(c(1:10),seq(from=12, to=20, by=2))
sft=pickSoftThreshold(input_mat_filt,powerVector = powers, verbose=5)
saveRDS(sft,file='wgcna/sft_power.RData')
sft<-readRDS("wgcna/sft_power.RData")

# plot diagnostics
png("wgcna/soft_threshold.png")
par(mfrow=c(1,2))
cex1=0.9
plot(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab='Soft Threshold (power)',ylab='Scale Free Topology Model Fit, Signed R^2',main=paste("Scale independence"))

text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, cex=cex1, col="red")
abline(h=0.90, col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main=paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")
dev.off()

# build the network
picked_power=10
temp_cor <- cor
cor <- WGCNA::cor
netwk <- blockwiseModules(input_mat_filt,power=picked_power,networkType = "signed", 
                          deepSplit=2,maxBlockSize = 20000, pamRespectsDendro = F, minModuleSize = 30, 
                          reassignThreshold = 0,mergeCutHeight = 0.25,
                          saveTOMs = T, saveTOMFileBase = "cult",numericLabels = T, verbose=3,corType = "bicor"
                          )
saveRDS(netwk,file='wgcna/netwk_bicor_blocksize.RData')
netwk<-readRDS("wgcna/netwk_bicor_blocksize.RData")

# look at results
mergedColors=labels2colors(netwk$colors)
png("wgcna/network_modules_bicor_blocksize.png")
plotDendroAndColors(netwk$dendrograms[[1]],mergedColors[netwk$blockGenes[[1]]], "Module colors", dendroLabels = FALSE,hang=0.03, addGuide = TRUE,guideHang = 0.05)
dev.off()

# print a text files of the genes/modules
module_df<-data.frame(gene_id=names(netwk$colors),colors=labels2colors(netwk$colors))
write_delim(module_df,file="wgcna/gene_modules.txt", delim='\t')

