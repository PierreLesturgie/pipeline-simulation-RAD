#!/bin/bash
#SBATCH -p workq
#SBATCH --time=6:10:00 #job time limit
#SBATCH --cpus-per-task=5 #ncpu on the same node



directory=`pwd`

cd $directory


module load system/Python-3.4.3


python3 Main_dig_parte_prima.py
	
