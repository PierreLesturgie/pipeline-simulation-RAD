#!/bin/bash
#SBATCH -p workq
#SBATCH --time=1:10:00 #job time limit
#SBATCH --cpus-per-task=10 #ncpu on the same node

module load system/Python-3.4.3

python3 Main_dig_parte_seconda.py
	


var1=$(wc -l posizioni.txt | cut -f1 -d' ')
rm *locus$var1.arp

