# Name: run_DGE_deseq_sunflower_inflo.R
# Author: EY
# Date: 01/24/2023 (edited 08/22/2023)
# Version:4.2.1
# Description: Will run the differential gene expression on the infloresence stages


library(dplyr)
library(DESeq2)
library(Glimma)
library(RUVSeq)
#devtools::install_github("zhangyuqing/sva-devel")
library(sva)

setwd('/home/ely67071/sunflower_inflo_dev_analysis/')

# read in the data matrix
summed_counts<-readRDS("/scratch/ely67071/sunflower_inflo_dev_data/collapsed_replicates_deseq_dataset.Rdata")
dim(summed_counts)



summed_counts$replicates


samples=c("10D_REP1_ATTACTCG", "20D_REP2_TCCGGAGA" ,"30D_REP2_CGCTCATT", "35D_REP1_GAGATTCC", 
          "HA_10D_2_ACCTTGGC", "HA_10D_3_ATATCTCG", "HA_20D_2_GCGCTCTA", 
          "HA_20D_3_AACAGGTT", "HA_30D_2_GGTGAACC", "HA_30D_3_CAACAATG", "HA_35D_2_TGGTGGCA", "HA_35D_3_AGGCAGAG")

dev_stage<-sub(".*([0-9]{2,2}D).*", "\\1",samples)
dev_stage <- as.factor(dev_stage)

metadata<-data.frame(samples, dev_stage)

# create the factors of interest
metadata$dev_stage<-factor(metadata$dev_stage)


# create the model
summed_counts<-DESeqDataSetFromMatrix(counts(summed_counts),colData = metadata, design=~0+dev_stage)

# pre-filter for reads where at least 3 samples have a count of 1 or higher
keep<-rowSums(counts(summed_counts)>=1)>=3
length(which(keep==1))
summed_counts_filt<-summed_counts[keep,]


write.csv(as.data.frame(counts(summed_counts)), file='deseq_results/raw_counts.csv')
as.data.frame(counts(summed_counts))


# run the DGE analysis
DESeq_dataset_results<-DESeq(summed_counts_filt,parallel=TRUE)

# save the R object
saveRDS(DESeq_dataset_results,file='DESeq_dataset_results_inflo.RData')

# load the dataset
#DESeq_dataset_results<-load("DESeq_dataset_results_inflo.RData")


# look at comparisons that can be made
resultsNames(DESeq_dataset_results)

# look at results for each treatment
results_10D<-results(DESeq_dataset_results, alpha=0.05,name='dev_stage10D')
results_20D<-results(DESeq_dataset_results, alpha=0.05,name='dev_stage20D')
results_30D<-results(DESeq_dataset_results, alpha=0.05,name='dev_stage30D')
results_35D<-results(DESeq_dataset_results, alpha=0.05,name='dev_stage35D')

write.csv(as.data.frame(results_10D), file='deseq_results/dev_stage10D.csv')
write.csv(as.data.frame(results_20D), file='deseq_results/dev_stage20D.csv')
write.csv(as.data.frame(results_30D), file='deseq_results/dev_stage30D.csv')
write.csv(as.data.frame(results_35D), file='deseq_results/dev_stage35D.csv')


result_10D_v_20D<-results(DESeq_dataset_results,contrast=c("dev_stage","10D","20D"),alpha=0.05,parallel=TRUE)
result_20D_v_30D<-results(DESeq_dataset_results,contrast=c("dev_stage","20D","30D"),alpha=0.05,parallel=TRUE)
result_30D_v_35D<-results(DESeq_dataset_results,contrast=c("dev_stage","30D","35D"),alpha=0.05,parallel=TRUE)
result_10D_v_35D<-results(DESeq_dataset_results,contrast=c("dev_stage","10D","35D"),alpha=0.05,parallel=TRUE)
result_10D_v_30D<-results(DESeq_dataset_results,contrast=c("dev_stage","10D","30D"),alpha=0.05,parallel=TRUE)
result_20D_v_35D<-results(DESeq_dataset_results,contrast=c("dev_stage","20D","35D"),alpha=0.05,parallel=TRUE)

write.csv(as.data.frame(result_10D_v_20D), file='deseq_results/pairwise/result_10D_v_20D.csv')
write.csv(as.data.frame(result_20D_v_30D), file='deseq_results/pairwise/result_20D_v_30D.csv')
write.csv(as.data.frame(result_30D_v_35D), file='deseq_results/pairwise/result_30D_v_35D.csv')


