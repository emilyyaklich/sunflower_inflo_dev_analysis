# Name: analyze DGE deseq sunflower inflo
# Author: EY (based off of code written by E. Dittmar)
# Date: 01/24/2023 (edited 08/22/2023)
# Version:4.2.1
# Description: Will analyze the output from the DESeq DGE for sunflower inflo stages
# need Functions.R written by ED


setwd('/home/ely67071/sunflower_inflo_dev_analysis/')

library(dplyr)
library(ggplot2)
library(UpSetR)
library(Glimma)
source("Functions.R")


# LOOK AT DATA BETWEEN TYPE
# read in the data
DEData_devstage<-ImportCSVs('deseq_results/',0.05)
# filter out significant results
mydataSig_devstage<-lapply(DEData_devstage, SigDEdf,PvaluesCol=7,CritP=0.05)


intersect(mydataSig_devstage$dev_stage10D[1], mydataSig_devstage$dev_stage20D[1])
mydataSig_devstage$dev_stage10D[1]





# see which genes overlap
SigOverlap_devstage<-GeneSets_four(mydataSig_devstage$dev_stage10D[1], mydataSig_devstage$dev_stage20D[1], mydataSig_devstage$dev_stage30D[1], mydataSig_devstage$dev_stage35D[1])
names(SigOverlap_devstage)
lapply(SigOverlap_devstage,function(x) {length(x$Gene)})


SigOverlapGraph_type<-lapply(mydataSig_type, function(x) {x$Gene})
png("cult/plots/type.png")
upset(fromList(SigOverlapGraph_type),order.by="freq",nsets=13,nintersects=20)
dev.off()





# LOOK AT DATA WHEN COMPARING SEQUENTIAL PAIRWISE
# read in the data
DEData_pairwise<-ImportCSVs('deseq_results/pairwise/',0.05)
# filter out significant results
mydataSig_pairwise<-lapply(DEData_pairwise,SigDEdf,PvaluesCol=7,CritP=0.05)
mydataSig_pairwise$result_10D_v_20D
mydataSig_pairwise$result_20D_v_30D
# see which genes overlap
SigOverlap_pairwise<-GeneSets(mydataSig_pairwise$result_10D_v_20D[1], mydataSig_pairwise$result_20D_v_30D[1],mydataSig_pairwise$result_30D_v_35D[1])
names(SigOverlap_pairwise)
lapply(SigOverlap_pairwise,function(x) {length(x$Gene)})

intersect(mydataSig_pairwise$result_10D_v_20D[1], mydataSig_pairwise$result_20D_v_30D[1])


#saveRDS(DEData_treat,SigOverlap_pairwise,file="Overlap.Rdata")

SigOverlapGraph_pairwise<-lapply(mydataSig_pairwise, function(x) {x$Gene})
lapply(mydataSig_pairwise, function(x) {x$Gene})
png("plots/sequential_pairwise_upset.png", res=215, width = 1200, height=1000)
upset(fromList(SigOverlapGraph_pairwise),order.by="freq",nsets=13,nintersects=20)
dev.off()





# w/ ruvseq


# LOOK AT DATA WHEN COMPARING SEQUENTIAL PAIRWISE
# read in the data
DEData_pairwise<-ImportCSVs('deseq_results/ruvseq/pairwise/',0.05)
# filter out significant results
mydataSig_pairwise<-lapply(DEData_pairwise,SigDEdf,PvaluesCol=7,CritP=0.05)
mydataSig_pairwise$result_10D_v_20D
mydataSig_pairwise$result_20D_v_30D
# see which genes overlap
SigOverlap_pairwise<-GeneSets(mydataSig_pairwise$result_10D_v_20D[1], mydataSig_pairwise$result_20D_v_30D[1],mydataSig_pairwise$result_30D_v_35D[1])
names(SigOverlap_pairwise)
lapply(SigOverlap_pairwise,function(x) {length(x$Gene)})

SigOverlapGraph_pairwise<-lapply(mydataSig_pairwise, function(x) {x$Gene})


png("plots/sequential_pairwise_upset_ruvseq.png", res=215, width = 1200, height=1000)
upset(fromList(SigOverlapGraph_pairwise),order.by="freq",nsets=13,nintersects=20)
dev.off()

















# with ruvS


# LOOK AT DATA WHEN COMPARING SEQUENTIAL PAIRWISE
# read in the data
DEData_pairwise<-ImportCSVs('deseq_results/ruvseq/ruvs/pairwise/',0.05)
# filter out significant results
mydataSig_pairwise<-lapply(DEData_pairwise,SigDEdf,PvaluesCol=7,CritP=0.05)
mydataSig_pairwise$result_10D_v_20D
mydataSig_pairwise$result_20D_v_30D
# see which genes overlap
SigOverlap_pairwise<-GeneSets(mydataSig_pairwise$result_10D_v_20D[1], mydataSig_pairwise$result_20D_v_30D[1],mydataSig_pairwise$result_30D_v_35D[1])
names(SigOverlap_pairwise)
lapply(SigOverlap_pairwise,function(x) {length(x$Gene)})

SigOverlapGraph_pairwise<-lapply(mydataSig_pairwise, function(x) {x$Gene})


png("plots/sequential_pairwise_upset_ruvseq.png", res=215, width = 1200, height=1000)
upset(fromList(SigOverlapGraph_pairwise),order.by="freq",nsets=13,nintersects=20)
dev.off()






# combatseq

# LOOK AT DATA WHEN COMPARING SEQUENTIAL PAIRWISE
# read in the data
DEData_pairwise<-ImportCSVs('deseq_results/combatseq/pairwise/',0.05)
# filter out significant results
mydataSig_pairwise<-lapply(DEData_pairwise,SigDEdf,PvaluesCol=7,CritP=0.05)
mydataSig_pairwise$result_10D_v_20D
mydataSig_pairwise$result_20D_v_30D
# see which genes overlap
SigOverlap_pairwise<-GeneSets(mydataSig_pairwise$result_10D_v_20D[1], mydataSig_pairwise$result_20D_v_30D[1],mydataSig_pairwise$result_30D_v_35D[1])
names(SigOverlap_pairwise)
lapply(SigOverlap_pairwise,function(x) {length(x$Gene)})

SigOverlapGraph_pairwise<-lapply(mydataSig_pairwise, function(x) {x$Gene})


png("plots/sequential_pairwise_upset_ruvseq.png", res=215, width = 1200, height=1000)
upset(fromList(SigOverlapGraph_pairwise),order.by="freq",nsets=13,nintersects=20)
dev.off()



