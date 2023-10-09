# Name: analyze DGE deseq sunflower inflo ruvseq
# Author: EY (based off of code written by E. Dittmar)
# Date: 08/29/2023
# Version:4.2.1
# Description: Will analyze the output from the DESeq DGE with ruvseq for sunflower inflo stages
# need Functions.R written by ED


setwd('/home/ely67071/sunflower_inflo_dev_analysis/')

library(dplyr)
library(ggplot2)
library(UpSetR)
library(Glimma)
source("Functions.R")


# now analyze 
# read in the data
DEData_pairwise_ruvs<-ImportCSVs('deseq_results/ruvseq/ruvs/pairwise/',0.05)
# filter out significant results
mydataSig_pairwise_ruvs<-lapply(DEData_pairwise_ruvs,SigDEdf,PvaluesCol=7,CritP=0.05)

# see which genes overlap
SigOverlap_pairwise_ruvs<-GeneSets(mydataSig_pairwise_ruvs$result_10D_v_20D[1], mydataSig_pairwise_ruvs$result_20D_v_30D[1],mydataSig_pairwise_ruvs$result_30D_v_35D[1])
names(SigOverlap_pairwise_ruvs)
lapply(SigOverlap_pairwise_ruvs,function(x) {length(x$Gene)})

SigOverlapGraph_pairwise_ruvs<-lapply(mydataSig_pairwise_ruvs, function(x) {x$Gene})

png("plots/sequential_pairwise_upset_ruvs.png", res=215, width = 1200, height=1000)
upset(fromList(SigOverlapGraph_pairwise_ruvs),order.by="freq",nsets=13,nintersects=20)
dev.off()



# now analyze 
# read in the data
DEData_pairwise_ruvr<-ImportCSVs('deseq_results/ruvseq/ruvr/pairwise/',0.05)
# filter out significant results
mydataSig_pairwise_ruvr<-lapply(DEData_pairwise_ruvr,SigDEdf,PvaluesCol=7,CritP=0.05)

# see which genes overlap
SigOverlap_pairwise_ruvr<-GeneSets(mydataSig_pairwise_ruvr$result_10D_v_20D[1], mydataSig_pairwise_ruvr$result_20D_v_30D[1],mydataSig_pairwise_ruvr$result_30D_v_35D[1])
names(SigOverlap_pairwise_ruvr)
lapply(SigOverlap_pairwise_ruvr,function(x) {length(x$Gene)})

SigOverlapGraph_pairwise_ruvr<-lapply(mydataSig_pairwise_ruvr, function(x) {x$Gene})

png("plots/sequential_pairwise_upset_ruvr.png", res=215, width = 1200, height=1000)
upset(fromList(SigOverlapGraph_pairwise_ruvr),order.by="freq",nsets=13,nintersects=20)
dev.off()




# now analyze RUVs data from LRT
DEData_pairwise_ruvs_lrt<-ImportCSVs('deseq_results/ruvseq/lrt/ruvs/',0.05)
# filter out significant results
mydataSig_pairwise_ruvs_lrt<-lapply(DEData_pairwise_ruvs_lrt,SigDEdf,PvaluesCol=6,CritP=0.05)

# see which genes overlap
SigOverlap_pairwise_ruvs_lrt<-GeneSets(mydataSig_pairwise_ruvs_lrt$result_10D_v_20D[1], mydataSig_pairwise_ruvs_lrt$result_20D_v_30D[1],mydataSig_pairwise_ruvs_lrt$result_30D_v_35D[1])
names(SigOverlap_pairwise_ruvs_lrt)
lapply(SigOverlap_pairwise_ruvs_lrt,function(x) {length(x$Gene)})

SigOverlapGraph_pairwise_ruvs_lrt<-lapply(mydataSig_pairwise_ruvs_lrt, function(x) {x$Gene})

png("plots/sequential_pairwise_upset_ruvs_lrt.png", res=215, width = 1200, height=1000)
upset(fromList(SigOverlapGraph_pairwise_ruvs_lrt),order.by="freq",nsets=13,nintersects=20)
dev.off()

