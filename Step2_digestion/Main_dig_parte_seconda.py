
# coding: utf-8

# In[ ]:

import library_digestion_in_silico
import read_seq


Liste_file = read_seq.lister_fichier_seq()


for i in Liste_file:
    library_digestion_in_silico.extract_single_digest(i, 'posizioni.txt', 100, 6)


