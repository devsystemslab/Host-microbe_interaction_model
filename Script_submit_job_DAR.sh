#!/bin/bash

for ct in Non_infection_specific_colonocyte Goblet Colon_SC Infection_specific_colonocyte_subtype_2 Infection_specific_colonocyte_subtype_1
do
    bsub -J DAR_${ct} -n 10 -o output-%J.out -e output-%J.err -R "rusage[mem=200G/host]affinity[core(1)*1:cpubind=core:membind=localprefer:distribute=any]" "/home/yuq22/miniconda3/envs/r-kernel/bin/Rscript /home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/ATAC_on_more_filtered_data/DAR/Script_identify_DAR_between_condition.R $ct"

done 