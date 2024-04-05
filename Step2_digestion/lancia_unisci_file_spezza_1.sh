#!/bin/bash
#SBATCH -p workq
#SBATCH --time=1:00:00 #job time limit
#SBATCH --cpus-per-task=1 #ncpu on the same node



directory=`pwd`


ini=0
fin=367

for c in `seq 1 250`;
do
mkdir "$directory"/set"$c"

cp lancia_unisci_file_spezza_2.sh "$directory"/set"$c"

cd "$directory"/set"$c"

sed -i 's|xxx|'$ini'|g' lancia_unisci_file_spezza_2.sh
sed -i 's|yyy|'$fin'|g' lancia_unisci_file_spezza_2.sh
sed -i 's|zzz|'$directory'|g' lancia_unisci_file_spezza_2.sh

sbatch lancia_unisci_file_spezza_2.sh

cd ..

ini=$(($ini+368))
fin=$(($fin+368))

done



