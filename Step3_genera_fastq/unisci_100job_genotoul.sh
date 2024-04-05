#!/bin/bash
#SBATCH -p workq
#SBATCH --time=4:10:00 #job time limit
#SBATCH --cpus-per-task=1 #ncpu on the same node


directory=`pwd`


cicli=250
ind=30



mkdir "$directory"/temp
for b in `seq 1 $cicli`;
do
for z in `seq 1 $ind`;
do
cp "$directory"/set"$b"/"$z".fasta "$directory"/temp/"$z".fasta."$b"
done
rm -r "$directory"/set"$b"
done

cd "$directory"/temp
for z in `seq 1 $ind`;
do
cat "$z".* > "$directory"/"$z".yy
done

cd ..

rm -r temp

for file in *.yy;
do
filename=$(basename -- "$file")
extension="${filename##*.}"
filename="${filename%.*}"
cp stampo_yy.sh aa."$file".sh

sed -i s/YY/"$file"/g  aa."$file".sh
sed -i s/QQ/"$filename"/g  aa."$file".sh
sbatch aa."$file".sh
done



