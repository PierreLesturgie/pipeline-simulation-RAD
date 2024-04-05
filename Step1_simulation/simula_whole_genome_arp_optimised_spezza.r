library(rlist)



crea_whole_genome_spezza<-function(percorso,genoma=1000000,versione="fsc26",segmenti=20){

nome<-sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(percorso))
path<-dirname(percorso)

##################################### PRENDO INPUT: LE SEQUENZE E LE POS DEGLI SNP
################### 
sequenze<-grep("_",readLines(files), value=T)
#versione=c("aaa")
#### estraggo solo la matrice dei dati da un file .arp: é la parte variabile delle future reads

n_ind<-length(sequenze)
if (versione=="fsc26") {temp_seq=strsplit(sequenze, split="\t")
} else {
temp_seq=strsplit(sequenze, split=" ")
}

matrice_nuc<-c()
for (i in 1:n_ind){
campo<-length(temp_seq[[i]])
matrice_nuc<-rbind(matrice_nuc,temp_seq[[i]][campo])
}
###############


##### prendo posizioni dal file arp per scrivere le parti conservate
pos_arl<-grep(",",readLines(files), value=T)
pos_arl_tmp=strsplit(pos_arl, split="#")  #### prende la linea snps e eventi ricombinazione
pos_arl_tmp_finali<-strsplit(pos_arl_tmp[[1]][2], split=", ") ### prendo quindi solo gli snps
n_snp_tot<-length(pos_arl_tmp_finali[[1]][])
pos_finali<-pos_arl_tmp_finali[[1]][]
pos_finali<-as.numeric(pos_finali)
#####
################################################### FINE LETTURA FILE: HO PRESO QUELLO DI CUI HO BISOGNO


###################################### ora creo tre liste con le cose di cui ho bisogno per scrive il genoma#
### ho bisogno di tre argomenti: pos_finali, matrice_nuc, genoma

############################ PRIMA LISTA: posizioni SNPs
##################### qui ceo pos_finali_list dove ho una list dove ogni elemento é un pos_finali. ci sono 'files_inali' elementi
n_seg<-segmenti ###ntervalli in cui dividere il numero tot di snp

l_seg<-floor(n_snp_tot/n_seg) ### lunghezza di tali intervalli

pos_finali_list<-list()

if (l_seg*n_seg<n_snp_tot) {
for (i in 1:n_seg){
pos_finali_list<-list.append(pos_finali_list,pos_finali[(1+((i-1)*l_seg)):(i*l_seg)])
}
pos_finali_list<-list.append(pos_finali_list,pos_finali[((n_seg*l_seg)+1):n_snp_tot])
files_finali<-n_seg+1
} else {
for (i in 1:n_seg){
pos_finali_list<-list.append(pos_finali_list,pos_finali[(1+((i-1)*l_seg)):(i*l_seg)])
}
files_finali<-n_seg
}


seq_finale_low<-seq(1,n_snp_tot,l_seg)
seq_finale_up<-c(seq_finale_low[2:length(seq_finale_low)]-1,n_snp_tot)

######################################## in seq_finali low eup ci sono gli estremi degli snp da prendere da matrice_nuc
##################################################################


######################################## SECONDA LISTA: SEQUENZE SNP
############################
matrice_nuc_list<-list()
qq<-strsplit(matrice_nuc, split="")
for (i in 1:files_finali) {
pezzo_temp<-c()
for (j in 1:n_ind){
pezzo_temp<-rbind(pezzo_temp, qq[[j]][seq_finale_low[i]:seq_finale_up[i]])
}
matrice_nuc_list<-list.append(matrice_nuc_list,pezzo_temp)
}
################################ qui in ogni elemento di matrice_nuc_list c'é la sequenza originaria di snps



######################################## TERZA LISTA: Parti da aggiungere conservate in ciascun dei file

intervalli=sort(c(1,genoma,pos_finali))

int_finale<-c()

for (i in 2:length(intervalli)) {
int_finale<-c(int_finale,intervalli[i]-intervalli[i-1])
}
int_finale<-int_finale-1
int_finale[1]<-int_finale[1]+1
int_finale[n_snp_tot+1]<-int_finale[n_snp_tot+1]+1
### int_finale contiene la dimensione di ciascuna parte. la loro somma + n_snp deve dare genoma Ha n_snp+1classi

###creo lista finale in cui ogni elemento contiene le lunghezze di ciasun frammento di parte conservata da aggiungere*

intervalli_list<-list()
for (i in 1:length(pos_finali_list)){
intervalli_list<-list.append(intervalli_list, int_finale[seq_finale_low[i]:seq_finale_up[i]])
}
intervalli_list[[length(pos_finali_list)]]<-c(intervalli_list[[length(pos_finali_list)]],int_finale[n_snp_tot+1])
########################################################################################################



####### ciclo finale in cui scrivo files_finali per individuo. il modo piu semplice é quindi cosi..

for (x in 1:files_finali){

n_snp<-length(pos_finali_list[[x]])
int_finale<-c(intervalli_list[[x]])
matrice_nuc<-matrice_nuc_list[[x]]

###### creo lista con le parti conservate. 


basi<-c("A","G","C","T")
finale<-list()

for (i in 1:(length(intervalli_list[[x]]))){

n_base<-int_finale[i] ###
parte_conservata<-c()
for (i in 1:n_base){
a=sample(basi,1)
parte_conservata<-c(parte_conservata,a)
}
finale<-list.append(finale,parte_conservata)
}

####################################### in finale ho tutte le regioni n_snp+1 conservate

#### scrivo nome dei file in una nuova lista
nome_output<-list()
for (i in 1:n_ind) {
nome_output<-list.append(nome_output,paste(path,"/", nome, "_", i , "_parte_", x, ".seq", sep=""))
}
#######

##### in quesyo ciclo scrivo un file alla volta, questo mi dovrebbe far risparmaire della memoria

for (j in 1:n_ind){

individui<-c()
for (i in 1:n_snp){
individui<-c(individui,finale[[i]],strsplit(matrice_nuc[j,], split="")[[i]])
}
#
if (x==files_finali) {individui<-c(individui,finale[[n_snp+1]])}

con <- file(nome_output[[j]], "w")
cat(individui, file=con, sep="")
cat("", file=con, sep="\n")
close(con)
rm(individui)
}

rm(finale)
rm(parte_conservata)
}

}



#setwd("C:/Users/smona/Documents/rad_simulati/whole_genome/new_optimised")
files <- list.files(pattern = "\\.arp$")
for (i in 1:length(files)) {

crea_whole_genome_spezza(files[i])

}

