#!/bin/bash
#SBATCH -p workq
#SBATCH --time=2:10:00 #job time limit
#SBATCH --cpus-per-task=1 #ncpu on the same node


directory=`pwd`

for b in `seq 1 250`;
do
cp * $directory/set"$b"
cd $directory/set"$b"
sbatch sim_crea_reads.sh
cd ..
done

