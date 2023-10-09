# Name: run_DGE_deseq_sunflower_inflo_ruvseq.R
# Author: EY
# Date: 08/29/2023
# Version:4.2.1
# Description: Will run the differential gene expression on the infloresence stages accounting for variation with RUVSeq


library(dplyr)
library(DESeq2)
library(Glimma)
library(RUVSeq)
library(ggplot2)
library(UpSetR)
library(edgeR)
source("Functions.R")


# read in and process data

setwd('/home/ely67071/sunflower_inflo_dev_analysis/')

# read in the data matrix
summed_counts<-readRDS("/scratch/ely67071/sunflower_inflo_dev_data/collapsed_replicates_deseq_dataset.Rdata")

samples=c("10D_REP1_ATTACTCG", "20D_REP2_TCCGGAGA" ,"30D_REP2_CGCTCATT", "35D_REP1_GAGATTCC", 
          "HA_10D_2_ACCTTGGC", "HA_10D_3_ATATCTCG", "HA_20D_2_GCGCTCTA", 
          "HA_20D_3_AACAGGTT", "HA_30D_2_GGTGAACC", "HA_30D_3_CAACAATG", "HA_35D_2_TGGTGGCA", "HA_35D_3_AGGCAGAG")

dev_stage<-sub(".*([0-9]{2,2}D).*", "\\1",samples)
dev_stage <- as.factor(dev_stage)

metadata<-data.frame(samples, dev_stage)

# create the factors of interest
metadata$dev_stage<-factor(metadata$dev_stage)



# pre-filter for reads where at least 3 samples have a count of 1 or higher
keep<-rowSums(counts(summed_counts)>=1)>=3
length(which(keep==1))
summed_counts_filt<-summed_counts[keep,]




# (1) RUVs (uses replicate samples to estimate factors of unwanted variation)

# get gene/transcript names
genes <- rownames(counts(summed_counts_filt))
# get a dataframe of the counts
count_matrix <- as.data.frame(counts(summed_counts_filt))
# create the RUVseq input matrix
input_set <- newSeqExpressionSet(as.matrix(count_matrix),phenoData = data.frame(dev_stage, row.names=colnames(count_matrix)))

# create a matrix which specifies the replicates
differences <- makeGroups(dev_stage)


# run RUVs
ruvs_output<-RUVs(input_set, genes, k=1, differences)
png("plots/pca_ruvs.png", res=215, width = 1200, height=1000)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plotPCA(ruvs_output,labels=FALSE, col=dev_stage)
legend("topright",inset=c(-0.4,0),legend=unique(dev_stage), fill=dev_stage)
dev.off()

# read into a DESeq dataset
summed_counts_filt_ruvs<-DESeqDataSetFromMatrix(countData = counts(ruvs_output),colData = pData(ruvs_output), design=~0+W_1+dev_stage)

summed_counts_filt_ruvs<-DESeqDataSetFromMatrix(countData = counts(ruvs_output),colData = pData(ruvs_output), design=~0+W_1+dev_stage)
pData(ruvs_output)

# run the DGE analysis
DESeq_dataset_results_ruvs<-DESeq(summed_counts_filt_ruvs,parallel=TRUE)

result_10D_v_20D_ruvs<-results(DESeq_dataset_results_ruvs,contrast=c("dev_stage","10D","20D"),alpha=0.05,parallel=TRUE)
result_20D_v_30D_ruvs<-results(DESeq_dataset_results_ruvs,contrast=c("dev_stage","20D","30D"),alpha=0.05,parallel=TRUE)
result_30D_v_35D_ruvs<-results(DESeq_dataset_results_ruvs,contrast=c("dev_stage","30D","35D"),alpha=0.05,parallel=TRUE)


write.csv(as.data.frame(result_10D_v_20D_ruvs), file='deseq_results/ruvseq/ruvs/pairwise/result_10D_v_20D.csv')
write.csv(as.data.frame(result_20D_v_30D_ruvs), file='deseq_results/ruvseq/ruvs/pairwise/result_20D_v_30D.csv')
write.csv(as.data.frame(result_30D_v_35D_ruvs), file='deseq_results/ruvseq/ruvs/pairwise/result_30D_v_35D.csv')



# now with RUVr (based on residuals from a first-pass GLM)


design <- model.matrix(~dev_stage, data=pData(input_set))
# first-pass GLM
y <- DGEList(counts=counts(input_set), group=dev_stage)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
res <- residuals(fit, type='deviance')

# run RUVr
ruvr_output <- RUVr(input_set,genes, k=1, res)

# plot a PCA of the output
png("plots/pca_ruvr.png", res=215, width = 1200, height=1000)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plotPCA(ruvr_output,labels=TRUE)
legend("topright",inset=c(-0.4,0),legend=unique(dev_stage), fill=dev_stage)
dev.off()
 
# create the DESeq object 
summed_counts_filt_ruvr<-DESeqDataSetFromMatrix(countData = counts(ruvr_output),colData = pData(ruvr_output), design=~0+W_1+dev_stage)

DESeq_dataset_results_ruvr<-DESeq(summed_counts_filt_ruvr,parallel=TRUE)


# look at comparisons that can be made
resultsNames(DESeq_dataset_results_ruvr)


result_10D_v_20D_ruvr<-results(DESeq_dataset_results_ruvr,contrast=c("dev_stage","10D","20D"),alpha=0.05,parallel=TRUE)
result_20D_v_30D_ruvr<-results(DESeq_dataset_results_ruvr,contrast=c("dev_stage","20D","30D"),alpha=0.05,parallel=TRUE)
result_30D_v_35D_ruvr<-results(DESeq_dataset_results_ruvr,contrast=c("dev_stage","30D","35D"),alpha=0.05,parallel=TRUE)


write.csv(as.data.frame(result_10D_v_20D_ruvr), file='deseq_results/ruvseq/ruvr/pairwise/result_10D_v_20D.csv')
write.csv(as.data.frame(result_20D_v_30D_ruvr), file='deseq_results/ruvseq/ruvr/pairwise/result_20D_v_30D.csv')
write.csv(as.data.frame(result_30D_v_35D_ruvr), file='deseq_results/ruvseq/ruvr/pairwise/result_30D_v_35D.csv')


# run RUVs data with LRT 

# normalize the counts
# set the design
design<-model.matrix(~0+W_1+dev_stage, data=pData(ruvs_output))
y<- DGEList(counts=counts(ruvs_output), group=dev_stage)
y<-calcNormFactors(y,method="TMM")


# estimate dispersion
y <-estimateGLMCommonDisp(y, design)
y <-estimateGLMTagwiseDisp(y,design)
plotBCV(y)

# fit the model
fit<-glmQLFit(y,design)

# use the output from the line below to design the contrasts 
colnames(fit)

res_10v20<-glmQLFTest(fit, contrast=c(0,-1,1,0,0))
summary(decideTests(res_10v20))
plotMD(res_10v20)
top_10v20<-topTags(res_10v20, n=Inf)
write.csv(as.data.frame(top_10v20), file='deseq_results/ruvseq/lrt/ruvs/result_10D_v_20D.csv')



res_20v30<-glmQLFTest(fit, contrast=c(0,0,-1,1,0))
summary(decideTests(res_20v30))
plotMD(res_20v30)
top_20v30<-topTags(res_20v30, n=Inf)
write.csv(as.data.frame(top_20v30), file='deseq_results/ruvseq/lrt/ruvs/result_20D_v_30D.csv')




res_30v35<-glmQLFTest(fit, contrast=c(0,0,0,-1,1))
summary(decideTests(res_30v35))
plotMD(res_30v35)
top_30v35<-topTags(res_30v35, n=Inf)
write.csv(as.data.frame(top_30v35), file='deseq_results/ruvseq/lrt/ruvs/result_30D_v_35D.csv')


