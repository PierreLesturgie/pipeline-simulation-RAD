library(rlist)



crea_reads_simcoal<-function(percorso, lunghezza_enzima=6, cov_mean=50, cov_sd=10, err_illumina=0.001, versione="fsc26") {  ###

nome<-sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(percorso))
path<-dirname(percorso)

sequenze<-grep("_",readLines(percorso), value=T)

#### estraggo solo la matrice dei dati da un file .arp: é la parte variabile delle future reads

n_ind<-length(sequenze)
if (versione=="fsc26") {
temp_seq=strsplit(sequenze, split="\t")
} else {
temp_seq=strsplit(sequenze, split=" ")
}

matrice_nuc<-c()
for (i in 1:n_ind){
campo<-length(temp_seq[[i]])
matrice_nuc<-rbind(matrice_nuc,temp_seq[[i]][campo])
}
reads_sim<- length(strsplit(matrice_nuc[1,],split="")[[1]])
###############



###### ora produco le reads per simulare coverage. per ogni cromosoma (e non individuo, quindi per ogni linea della matrice_nuc)
###### replico delle reads, il cui numero sarà estratto da una normale. Se voglio simulare un coverage medio di 40X, ogni allele
###### avra un coverage di 20X. Metto una sd di 10, tutti e due parametri da mettere nella funzione eventuale
###### creo una matrice per individuo, prendendo a caso due linee di matrice_nuc per replicarle
###### conviene quindi scrivere per ogni locus n_ind nuove matrici, perché le devo salvare in file diversi ed é + pratico
###### simulo almeno 1 reads per locus, per una questione di semplicità con gli indici

mean_cov<-cov_mean
sd_cov<-cov_sd



### creo locus con tutte le reads
locus<-c()
indici_alleli<-c() ### contiene quante reads ci sono per individuo
for (i in 1:nrow(matrice_nuc)){
aggiungi<-floor(abs(rnorm(1,mean_cov,sd_cov)))
if (aggiungi==0){aggiungi=1}
#print(aggiungi)
indici_alleli<-c(indici_alleli, aggiungi)
for (j in 1:aggiungi){
locus<-rbind(locus,matrice_nuc[i,])
}
}

######


### mi scrivo gli indici per sapere a che cromosoma corrisponde ciascuna reads
born_inf<-c()
born_sup<-c()
for (z in 1:n_ind){
####ciclo per prendere i limiti di ciascuna pop dove vedere se lo SNP c'é in almeno un individuo
if (z==1) {
born_inf=cbind(born_inf,1)
born_sup=cbind(born_sup,indici_alleli[z])
}
else {

born_inf=cbind(born_inf,(sum(indici_alleli[1:(z-1)])+1))
born_sup=cbind(born_sup,(sum(indici_alleli[1:z])))
}
}

####


#### creo la lista dove ogni elemento é un individo (i due alleli son stati creati indipendentemente)
finale<-list()
for (z in 1:(n_ind/2)){
mat_temp<-locus[born_inf[(z*2-1)]:born_sup[(z*2)],]
finale<-list.append(finale,mat_temp)
}
##### finale é una lista con n_ind/2 individui (n_ind é il numero di cromosomi).

###########
###########
###########

########################################################## aggiungo parte conservata a ciascun locus +  mutazioni
reads_totali_per_ind<-list()

for (m in 1:length(finale)){



reads_intera<-c()
for (i in 1:length(finale[[m]])){
reads_intera<-rbind(reads_intera,strsplit(finale[[m]][i], split="")[[1]])
}

locus_finale<-reads_intera ##### creo data frame reads

########## aggiungo mutazioni randon, simulando errore Illumina

epsilon=err_illumina #### il valore é espresso come probabilità per base. posso metterlo come parametro dell'eventuale funzione
n_tot_mut<-floor(dim(locus_finale)[1]*dim(locus_finale)[2]*epsilon)
if (n_tot_mut==0) {n_tot_mut=1}

for (z in 1:n_tot_mut){

i<-floor(runif(1,1,dim(locus_finale)[1]))
j<-floor(runif(1,1,dim(locus_finale)[2]))

basi<-c("A","G","C","T")
a=sample(basi,1)
while (a == locus_finale[i,j]) {
a=sample(basi,1)
}

locus_finale[i,j]=a
}


reads_totali_per_ind<-list.append(reads_totali_per_ind,locus_finale)
}

########  reads_totali_per_ind  é una lista dove ogni elemento sono le reads con mut per individuo da scrivere in file fastq


####### scrivi ogni elemento della lista in un file diverso come fosse una vera read in fastq

#### scrivo nome dei file in una nuova lista
nome_output<-list()
for (i in 1:(n_ind/2)) {
nome_output<-list.append(nome_output,paste(path,"/", nome, "_", i , ".fq", sep=""))
}

for (i in 1:(n_ind/2)) {
con <- file(nome_output[[i]], "w")

for (j in 1:nrow(reads_totali_per_ind[[i]])){
cat("@", file=con)
cat(rnorm(1,12,13), file=con)
cat(rnorm(1,12,13), file=con)
cat("/1", file=con, sep="\n")
cat(reads_totali_per_ind[[i]][j,((lunghezza_enzima+1):reads_sim)], file=con, sep="")
cat("", file=con, sep="\n")
cat("+", file=con, sep="\n")
cat(rep("F",(reads_sim-lunghezza_enzima-10)), file=con, sep="")
cat(rep("G",10), file=con, sep="")
cat("", file=con, sep="\n")
}
close(con)
}
}


files <- list.files(pattern = "\\.fine$")
for (i in 1:length(files)) {

crea_reads_simcoal(files[i])

}
