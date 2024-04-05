#!/bin/bash
#SBATCH -p workq
#SBATCH --time=40:10:00 #job time limit
#SBATCH --cpus-per-task=1 #ncpu on the same node

module load system/R-3.5.1


directory=`pwd`
ind=30

Rscript --vanilla simula_reads_funzione_da_wholegenome.r


for z in `seq 1 $ind`;
do
cat *_"$z".fq > "$directory"/"$z".fasta 
done
cd ..
