#!/bin/bash
#SBATCH --job-name=prepare_DEseq_inputs
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=1gb
#SBATCH --time=40:00:00
#SBATCH --output=ErrorFiles/prepare_DEseq_inputs.%j.out
#SBATCH --error=ErrorFiles/prepare_DEseq_inputs.%j.error
#SBATCH --mail-user=ely67071@uga.edu
#SBATCH --mail-type=ALL



for fn in /scratch/ely67071/sunflower_inflo_dev_data/gene_count_data/*.tab; do
    awk -F '\t' 'NR>4 { print $1, $4 }' "$fn" >tmp_file 
    mv tmp_file "$fn"
done    
