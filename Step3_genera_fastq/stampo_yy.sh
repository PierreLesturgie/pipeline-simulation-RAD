#!/bin/bash
#SBATCH -p workq
#SBATCH --time=10:10:00 #job time limit
#SBATCH --cpus-per-task=1 #ncpu on the same node


#
directory=`pwd`


awk '{printf("%s%s",$0,(NR%4==0)?"\n":"\0")}' YY | sort -R | tr "\0" "\n" > QQ.fq

rm YY


