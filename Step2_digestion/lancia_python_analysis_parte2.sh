#!/bin/bash
#SBATCH -p workq
#SBATCH --time=1:10:00 #job time limit
#SBATCH --cpus-per-task=1 #ncpu on the same node


directory=`pwd`

ind=60
for b in `seq 1 $ind`;
do
mkdir "$directory"/seq_"$b"
cp "$directory"/"$b".fasta  "$directory"/seq_"$b"
cp "$directory"/*.py "$directory"/seq_"$b"
cp "$directory"/*.sh "$directory"/seq_"$b"
cp "$directory"/posizioni.txt "$directory"/seq_"$b"
cd "$directory"/seq_"$b"
sed -i 's|xxxx|seq_'$b'|g' read_seq.py
sbatch Python_Analysis_parte_seconda.sh
cd ..
done



