# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 15:16:57 2025

@author: sinan.oz
"""

import numpy as np
import networkx as nx
import pickle
import pandas as pd 
import math
import matplotlib.pyplot as plt
from operator import itemgetter

data1=pd.read_excel('adjacency_matrix_dodecanes.xlsx')
data2=pd.read_excel('atom_pair_dodecanes.xlsx')
data3=pd.read_excel('topological_torsion_dodecanes.xlsx')
data4=pd.read_excel('Avalon_dodecanes.xlsx')
data5=pd.read_excel('MACCS_keys_dodecanes.xlsx')
data6=pd.read_excel('Morgan_circular_dodecanes.xlsx')


#Tanimoto coefficient values versus the number of dodecanes.
adjacency_matrix_dodecanes=data1[0].tolist()
ordered_adjacency_matrix_dodecanes=sorted(adjacency_matrix_dodecanes)
atom_pair_dodecanes=data2[0].tolist()
ordered_atom_pair_dodecanes=sorted(atom_pair_dodecanes)
topological_torsion_dodecanes=data3[0].tolist()
ordered_topological_torsion_dodecanes=sorted(topological_torsion_dodecanes)
Avalon_dodecanes=data4[0].tolist()
ordered_Avalon_dodecanes=sorted(Avalon_dodecanes)
MACCS_keys_dodecanes=data5[0].tolist()
ordered_MACCS_keys_dodecanes=sorted(MACCS_keys_dodecanes)
Morgan_circular_dodecanes=data6[0].tolist()
ordered_Morgan_circular_dodecanes=sorted(Morgan_circular_dodecanes)


number=[]
for i in range (0,62835):
    number.append(i)
    
# Plotting
ax = plt.subplot(111)
plt.plot(number, ordered_adjacency_matrix_dodecanes, label="Based on adjacency matrix ") 
plt.plot(number, ordered_atom_pair_dodecanes, label="Based on atom pair")
plt.plot(number, ordered_topological_torsion_dodecanes, label="Based on topological torsion") 
plt.plot(number, ordered_Avalon_dodecanes, label="Based on Avalon") 
plt.plot(number, ordered_MACCS_keys_dodecanes, label="Based on MACCS keys") 
plt.plot(number, ordered_Morgan_circular_dodecanes, label="Based on Morgan circular") 
ax.legend(bbox_to_anchor=(1.1, 1.05)) 
plt.xlabel('Number of graphs to be compared')
plt.ylabel('Tanimoto coefficient')
plt.suptitle('For Dodecanes')
#

