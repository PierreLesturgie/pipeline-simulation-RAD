
# coding: utf-8

# In[ ]:

import os
import pandas as pan
import re 
import sys
import csv

csv.field_size_limit(sys.maxsize)

def single_digest(filename, motif, file_output):

    with open(filename, newline='') as inputfile:
        results = list(csv.reader(inputfile))

    DNA_test = open(filename).read()
    motif=str(motif)
    aaa=[]
    for match in re.finditer(motif, DNA_test):
        aaa.append(match.span())
        #print(match.span(), match.group())

    l=len(aaa)
    print(l)
    ris=[]
    for i in range(l):
        ris.append(aaa[i][0])

    with open(file_output, 'w') as f:
        f.write('\n'.join(str(i) for i in ris))
        f.close()
    return(l)




def extract_single_digest(sequences, positions,lunghezza, lungh_enzima):
    
    with open(sequences, newline='') as inputfile: ### qui apro il file contenente le sequenze consensus. Ne devo avere una per linea
        results = list(csv.reader(inputfile))

    DNA_test = open(sequences).read().split() ### le trasformo in matrice di dati

    with open(positions, newline='') as input2: ### apro le posizioni da estrarre nel multiallineamento
        res2= list(csv.reader(input2))
    posizioni = pan.DataFrame(res2) ### le trasformo in un data frame

    loci=len(posizioni) ### questo é quindi il numero totale di loci
    print(loci)
    l = len(DNA_test)  #questo é il numero di individui
    print(l)
    l_frammento=lunghezza  ### definisco la lunghezza del frammento da estrarre
    l_enz=lungh_enzima

    loci_fin={} ### creo il mio dizionario. Ogni elemento sarà un locus estratto dalle mie sequenze
    for j in range(loci):
        DNA_locus=[]
        for i in range(l):    
            s= DNA_test[i]                      
            inizio=int(posizioni[0][j]) ## definisco il punto di inzio e fine da estrarre
            fine=int(posizioni[0][j])+l_frammento+l_enz
            restriction_site= s[inizio:fine]
            DNA_locus.append(restriction_site) ### creo il locus j sui miei L individui
        loci_fin[j]=DNA_locus ### metto il locus j nel mio dizionario


    nome_out,estensione=os.path.splitext(sequences)
    for j in range(loci): ### vado a scrivere un files per locus
        filename=nome_out+"_locus"+str(j)+".arp"
        f = open(filename, "w", newline='')
        print(filename)
        for i in range(l): ### in ogni file ci metto tutte le entries del mio dizionario
            f.write("1_1\t"+"1\t") ### scrivo il file stile arelquin per farlo leggere agli altri script
            f.write(loci_fin[j][i])
            f.write("\n")
        f.close()
    return(0)




def double_digest(filename, motif1, motif2, taglia_low, taglia_up, file_output):

    with open(filename, newline='') as inputfile:
        results = list(csv.reader(inputfile))

    DNA_test = open(filename).read()
    enz1=[]
    for match in re.finditer(motif1, DNA_test): ### con questo loop trovo tutte le occorrenze del mio pattern
        enz1.append(match.span()) ### metto solo le posizioni di inizio e fine del mio pattern in aaa
        print(match.span(), match.group())

    l1=len(enz1) ### mi dice quante volte ho trovato il pattern
    print(l1)
    ris_enz1=[]
    for i in range(l1): ### estraggo solo la posizione di inizio del mio pattern
        ris_enz1.append(enz1[i][0])
    print(ris_enz1)


    enz2=[]
    for match in re.finditer(motif2, DNA_test): ### con questo loop trovo tutte le occorrenze del mio pattern
        enz2.append(match.span()) ### metto solo le posizioni di inizio e fine del mio pattern in aaa
        print(match.span(), match.group())

    l2=len(enz2) ### mi dice quante volte ho trovato il pattern
    print(l2)
    ris_enz2=[]
    for i in range(l2): ### estraggo solo la posizione di inizio del mio pattern
        ris_enz2.append(enz2[i][0])
    print(ris_enz2)

    print(len(ris_enz2))

    frammenti_low=[]
    frammenti_up=[]
    enzima1=[] ### qui metto l'enzima che taglia a valle. Mi serve per estrarre quel frammento di DNA che poi incollo a monte
    enzima2=[] ### qui metto l'enzima che taglia a monte. mi serve per capire quanti nt estrarre dopo (affinché sia uguale la lungghezza indipendentemente dalla lunghezza dei due siti)
    for i in range(len(ris_enz1)):
        for j in range(len(ris_enz2)):
            if (int(ris_enz1[i]) - int(ris_enz2[j])) > taglia_low and (int(ris_enz1[i]) - int(ris_enz2[j])) < taglia_up:
                #print(int(ris_enz1[i]))
                #print(int(ris_enz2[j]))
                frammenti_up.append(ris_enz1[i])
                frammenti_low.append(ris_enz2[j])
                enzima1.append(motif1) ### metto quello a valle (in questo calo l'1). E' la sua lunghezza che mi dice quante basi dovro poi estrarre
                enzima2.append(motif2)
    for i in range(len(ris_enz2)):
        for j in range(len(ris_enz1)):
            if (int(ris_enz2[i]) - int(ris_enz1[j])) > taglia_low and (int(ris_enz2[i]) - int(ris_enz1[j])) < taglia_up:
                #print(int(ris_enz1[j]))
                #print(int(ris_enz2[i]))
                frammenti_up.append(ris_enz2[i])
                frammenti_low.append(ris_enz1[j])
                enzima1.append(motif2) ### metto quello a valle (in questo calo il 2). E' la sua lunghezza che mi dice quante basi dovro poi estrarre
                enzima2.append(motif1)

    with open(file_output, 'w') as f:
        for i in range(len(frammenti_low)):
            f.write((str(frammenti_low[i])))
            f.write("\t")
            f.write((str(frammenti_up[i])))
            f.write("\t")
            f.write((str(enzima1[i])))
            f.write("\t")
            f.write((str(enzima2[i])))
            f.write("\n")
        f.close()
        
        
        
        
def extract_double_digest(sequences, positions, lunghezza):
    
    with open(sequences, newline='') as inputfile: ### qui apro il file contenente le sequenze consensus. Ne devo avere una per linea
        results = list(csv.reader(inputfile))

    DNA_test = open(sequences).read().split() ### le trasformo in matrice di dati

    with open(positions, newline='') as input2: ### apro le posizioni da estrarre nel multiallineamento
        res2= list(list(csv.reader(input2, delimiter='\t')))
    posizioni = pan.DataFrame(res2) ### le trasformo in un data frame


    loci=len(posizioni) ### questo é quindi il numero totale di loci
    print(loci)
    l = len(DNA_test)  #questo é il numero di individui
    print(l)
    l_frammento=lunghezza  ### definisco la lunghezza del frammento da estrarre


    loci_fin={} ### creo il mio dizionario. Ogni elemento sarà un locus estratto dalle mie sequenze
    for j in range(loci): ### con questo for prendo il frammento a valle del primo sito di taglio
        DNA_locus=[]
        for i in range(l):    
            s= DNA_test[i]                      
            inizio=int(posizioni[0][j]) ## definisco il punto di inzio e fine da estrarre
            fine=int(posizioni[0][j])+l_frammento+len(posizioni[3][j])
            restriction_site= s[inizio:fine]
            DNA_locus.append(restriction_site) ### creo il locus j sui miei L individui
        loci_fin[j]=DNA_locus ### metto il locus j nel mio dizionario

    
    loci_fin_valle={}    
    for j in range(loci): ### con questo for prendo il frammento a valle del secondo sito di taglio:
        DNA_locus_valle=[]  ## solo la lunghezza del sito di restrizione corrispondente, che ho nel file di input
        for i in range(l):    
            s= DNA_test[i]                      
            inizio=int(posizioni[1][j]) ## definisco il punto di inzio e fine da estrarre
            fine=int(posizioni[1][j])+len(posizioni[2][j]) ## la fine si base su lunghezza sito restrizione a valle della digestione (che ho szcritto nel file di input)
            restriction_site_valle= s[inizio:fine]
            DNA_locus_valle.append(restriction_site_valle) ### creo il locus j sui miei L individui
        loci_fin_valle[j]=DNA_locus_valle ### metto il locus j nel mio dizionario

    nome_out,estensione=os.path.splitext(sequences)
    for j in range(loci): ### vado a scrivere un files per locus
        filename=nome_out+"_locus"+str(j)+".arp"
        f = open(filename, "w", newline='')
        print(filename)
        for i in range(l): ### in ogni file ci metto tutte le entries del mio dizionario
            f.write("1_1\t"+"1\t") ### scrivo il file stile arelquin per farlo leggere agli altri script
            f.write(loci_fin_valle[j][i]) ### per comodita di calcolo (dovuta agli script calcolo sumstat+sfs), attacco il sito di restrizione valle a monte della sequenza!!!!cosi nei file che scrivo ho i due siti in continuità e posso utilizzare script fatti per le sim da fastsimcoal
            f.write(loci_fin[j][i])
            f.write("\n")
        f.close()

    return(0)

