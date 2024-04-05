#!/bin/bash
#SBATCH -p workq
#SBATCH --time=1:10:00 #job time limit
#SBATCH --cpus-per-task=1 #ncpu on the same node


directory=`pwd`

for b in `seq 1 5`;
do
mkdir "$directory"/set"$b"
cp * $directory/set"$b"
cd $directory/set"$b"
sbatch sim_whole_reads.sh
cd ..
done

