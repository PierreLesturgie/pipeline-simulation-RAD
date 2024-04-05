# Pipeline for simulating RADseq loci 

## The pipeline works in **3 steps** : 



### - Step 1. Simulation of N haploid individuals from a demographic model 
Folder ***Step1_simulation***

Five files are needed, they should be present in the same folder. The demographic model is specified in the file *.par*. 

#### To perform the simulations: 
1. **Execute** *lancia_6_job.sh*. This splits the job in different nodes of a computer cluster. Each job corresponds to an independent chromosome of the size specified in the *.par* file).
2. After the simulations has finished and that the full sequence has been written, **execute** *crea_fasta_finali.sh*. This merges all files to create the final *N* files, each containing one haploid sequence (*.fasta*).

**NOTE:** in the *simula_whole_genome_arp_optimised_spezza.r* you need to specify the length of the sequence to be obtained (which needs to correspond to the *.par* file).


### - Step 2. Digestion 
Folder ***Step2_digestion***

To perform digestion: 
1. **Execute** *Python_Analysis_new_parte_prima.sh*. This extracts the cut position from one sequence specified in the file *Main_dig_parte_prima.py*. In this script, there is only one function with three arguments: **(1) name of the sequence to analyze**; **(2) motif of the restriction enzyme**, and **(3) name of the output file storing the coordinates of the cutting points**.
2. **Execute** *lancia_python_analysis_parte2.sh*. This creates N folder containing the loci for each individual. The number of base pairs extracted after the cut points is specified in the file *Main_dig_parte_seconda.py* (3rd arguments). The 4th argument is the length of the restriction enzyme motif.
3. **Execute** *lancia_unisci_file_spezza_1.sh* to concatenate the same locus for all individuals in a single file. To do this, the sript specifies **(1) the number of folders to create (called *set*)**, here 250; **(2) how many loci put in each folder**, here 367 (which is done by executing *lancia_unisci_file_spezza_2.sh*). In each folder there will be 367 files called locus*.fine. Each file contains N lines corresponding to the sequence of the * locus in the N individuals.

**NOTE:** You need to specify where are the sequences generated in step 1 in the file read_seq.py (repertoire =).
  
  
### - Step 3. 
files in the folder script_genera_fastq_brand_new. Files need to be executed where the set* folders are. First, execute lancia_100job.sh which will run the R script simula_reads_funzione_da_wholegenome.r to all loci. In the R function within the script (crea_reads_simcoal), four arguments can be set: 1) the emplacement of the file (which need not to be changed), the length of the restriction enzyme used, the mean and standard deviation of the normal distribution from which to extract the N° of reads per individual per locus, the Illumina error to add to each reads, the version of fastsimcoal used for simulations. This function will create a fastq file for each diploid individual (merging the same couple of N haploid individuals for each locus form the *.fine files) containing the reads for all locus in the set* folder. Finally, execute the file unisci_100job_genotoul.sh (specifying the N° of diploid individuals and the number of set* folders) which will merge all loci for each individuals and randomly shuffle the reads to mimic a real output of an Illumina machine. The final result is one fastq files containing all loci for each diploid individual. These files will have the depth of coverage specified in the R function crea_reads_simcoal and can be subsampled to lower it (since the reads have been shuffled). 
