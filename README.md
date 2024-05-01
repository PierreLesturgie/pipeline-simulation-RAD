# Pipeline for simulating RADseq loci 
----

#### Author: Stefano Mona (stefano.mona@mnhn.fr), Ecole Pratique des Hautes Etudes - Paris Sciences Lettres

This is a pipeline written for: 

Mona, S., Benazzo, A., Delrieu-Trottin, E., & Lesturgie, P. (2023). Population genetics using low coverage RADseq data in non-model organisms: biases and solutions. https://doi.org/10.22541/au.168252801.19878064/v1

Please cite this paper when using the pipeline. 

----


The pipeline works in **3 steps** : **(1)** Simulation; **(2)** Digestion; **(3)** Generation of FASTQ files. 


## (1) Simulation
Folder ***Step1_simulation***

Simulation of N haploid individuals from a demographic model 

Five files are needed, they should be present in the same folder. The demographic model is specified in the file *.par*. 

#### To perform the simulations: 
1. **Execute** *lancia_6_job.sh*. This splits the job in different nodes of a computer cluster. Each job corresponds to an independent chromosome of the size specified in the *.par* file.
2. After the simulations has finished and that the full sequence has been written, **execute** *crea_fasta_finali.sh*. This merges all files to create the final *N* files, each containing one haploid sequence (*.fasta*).

**NOTE:** in the *simula_whole_genome_arp_optimised_spezza.r* you need to specify the length of the sequence to be obtained (which needs to correspond to the *.par* file).


## (2) Digestion 
Folder ***Step2_digestion***

To perform digestion: 
1. **Execute** *Python_Analysis_new_parte_prima.sh*. This extracts the cut position from one sequence specified in the file *Main_dig_parte_prima.py*. In this script, there is only one function with three arguments: **(1) name of the sequence to analyze**; **(2) motif of the restriction enzyme**, and **(3) name of the output file storing the coordinates of the cutting points**.
2. **Execute** *lancia_python_analysis_parte2.sh*. This creates N folder containing the loci for each individual. The number of base pairs extracted after the cut points is specified in the file *Main_dig_parte_seconda.py* (3rd arguments). The 4th argument is the length of the restriction enzyme motif.
3. **Execute** *lancia_unisci_file_spezza_1.sh* to concatenate the same locus for all individuals in a single file. To do this, the sript specifies **(1) the number of folders to create (called *set*)**, here 250; **(2) how many loci put in each folder**, here 367 (which is done by executing *lancia_unisci_file_spezza_2.sh*). In each folder there will be 367 files called locus*.fine. Each file contains N lines corresponding to the sequence of the * locus in the N individuals.

**NOTE:** You need to specify where are the sequences generated in step 1 in the file read_seq.py (repertoire =).
  
  
## (3) Generate fastq

Folder ***Step3_genera_fastq***
To Generate the fastq: 
1. **Execute** *lancia_100job.sh*. This runs the R script *simula_reads_funzione_da_wholegenome.r* for all loci. In the R function within the script (crea_reads_simcoal), six arguments can be set: **(1)** the emplacement of the file [*percorso*] (which need not to be changed); **(2)** the length of the restriction enzyme used [*lunghezza_enzima*]; **(3)** the mean and **(4)** standard deviation of the normal distribution from which to extract the N° of reads per individual per locus [*cov_mean*] & [*cov_sd*]; **(5)** the Illumina error to add to each reads [*err_illumina*] and **(6)** the version of fastsimcoal used for simulations [*versione*]. This function will create a fastq file for each diploid individual (merging the same couple of N haploid individuals for each locus form the *.fine* files) containing the reads for all locus in the *set* folder.
2. **Execute** *unisci_100job_genotoul.sh* (specifying the N° of diploid individuals and the number of set* folders). This merges all loci for each individuals and randomly shuffle the reads to mimic a real output of an Illumina machine.

The final result is one fastq files containing all loci for each diploid individual. These files will have the depth of coverage specified in the R function crea_reads_simcoal and can be subsampled to lower it (since the reads have been shuffled). 

 **NOTE:** Scripts need to be executed where the *set* folders are.
