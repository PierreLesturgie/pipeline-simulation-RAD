
# coding: utf-8

# In[4]:


import os
import pandas as pan
from os import listdir


repertoire = "/home/vnicolas/work/ste_prova_radsim/for_digestion/xxxx/"  # donne le main directory ou se trouve le script Dropout.py (ce script)


def lister_fichier_seq ():
    liste_fichier_seq = []
    for nom_fichier in listdir(repertoire):
        if nom_fichier.endswith(".fasta"):
            liste_fichier_seq.append(repertoire +"/" +  nom_fichier )
            
    return liste_fichier_seq

