#!/bin/bash
#SBATCH --job-name=load_GC_data_and_sum_reps
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=100gb
#SBATCH --time=40:00:00
#SBATCH --output=hpc_output/load_GC_data_and_sum_reps.%j.out
#SBATCH --error=hpc_output/load_GC_data_and_sum_reps.%j.error
#SBATCH --mail-user=ely67071@uga.edu
#SBATCH --mail-type=ALL


# load R
ml R/4.2.1-foss-2020b

Rscript load_GC_data_and_sum_reps.R

#Rscript run_edgeR.R

#Rscript load_GC_and_sum_reps_deseq.R
