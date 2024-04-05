#!/bin/bash
#SBATCH -p workq
#SBATCH --time=40:00:00 #job time limit
#SBATCH --cpus-per-task=1 #ncpu on the same node

module load system/R-3.5.1
module load bioinfo/fsc2603


directory=`pwd`
filename=$(basename -- "*.par")
cartella=$(basename $filename .par)


fsc26 -i $filename  -n1


cp simula_whole_genome_arp_optimised_spezza.r "$cartella"
cd "$cartella"
Rscript --vanilla simula_whole_genome_arp_optimised_spezza.r

