#!/bin/bash
#SBATCH -p workq
#SBATCH --time=1:00:00 #job time limit
#SBATCH --cpus-per-task=10 #ncpu on the same node

directory=`pwd`



for d in `seq 1 60`;
do
for e in `seq xxx yyy`;
do
mv zzz/seq_"$d"/*locus$e.arp "$directory"
done
done


for f in `seq xxx yyy`; 
do
cat *locus$f.arp > locus$f.fine
rm *locus$f.arp
done

