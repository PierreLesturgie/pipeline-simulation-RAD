#!/bin/bash
#SBATCH -p workq
#SBATCH --time=10:00:00 #job time limit
#SBATCH --cpus-per-task=5 #ncpu on the same node


ind=60
cicli=5
directory=`pwd`


filemod=$(basename -- "*.par")
modello=$(basename $filemod .par)


for j in `seq 1 $cicli`;
do
cd set"$j"/"$modello"

for b in `seq 1 $ind`;
do
filename=$(basename -- "*.arp")
cartella=$(basename $filename .arp)
cat "$cartella"_"$b"_*.seq > "$b".fine
cat "$b".fine | tr -d '\n' > "$b".fasta
done
mv *.fasta "$directory"/set"$j"
cd "$directory"

done

mkdir temp
for j in `seq 1 $cicli`;
do

for b in `seq 1 $ind`;
do
mv "$directory"/set"$j"/"$b".fasta "$directory"/temp/"$b"_"$j".bene
done

done

for b in `seq 1 $ind`;
do
cat "$directory"/temp/"$b"_*.bene | tr -d '\n' > "$directory"/"$b".fasta
done



