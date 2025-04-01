#!/bin/bash

for ct in Colonocyte_4 Colonocyte_2 Colonocyte_1 Colonocyte_3 Goblet Colon_SC
do
    bsub -J GRN_${ct} -n 10 -o output-%J.out -e output-%J.err -R "rusage[mem=200G/host]affinity[core(1)*1:cpubind=core:membind=localprefer:distribute=any]" "/home/yuq22/miniconda3/envs/r-kernel/bin/Rscript /home/yuq22/Bacteria_TRM/salmonella/multiome_analysis/GRN_on_more_filtered_data/GRN_per_cell_type/Script_build_GRN.R $ct"

done 
