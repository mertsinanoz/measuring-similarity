# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 10:54:00 2025

@author: sinan.oz
"""

import numpy as np
import networkx as nx
import pickle
import pandas as pd 
import math
import matplotlib.pyplot as plt
from operator import itemgetter

data1=pd.read_excel('adjacency_matrix_undecanes.xlsx')
data2=pd.read_excel('atom_pair_undecanes.xlsx')
data3=pd.read_excel('topological_torsion_undecanes.xlsx')
data4=pd.read_excel('Avalon_undecanes.xlsx')
data5=pd.read_excel('MACCS_keys_undecanes.xlsx')
data6=pd.read_excel('Morgan_circular_undecanes.xlsx')


#Tanimoto coefficient values versus the number of undecanes.
adjacency_matrix_undecanes=data1[0].tolist()
ordered_adjacency_matrix_undecanes=sorted(adjacency_matrix_undecanes)
atom_pair_undecanes=data2[0].tolist()
ordered_atom_pair_undecanes=sorted(atom_pair_undecanes)
topological_torsion_undecanes=data3[0].tolist()
ordered_topological_torsion_undecanes=sorted(topological_torsion_undecanes)
Avalon_undecanes=data4[0].tolist()
ordered_Avalon_undecanes=sorted(Avalon_undecanes)
MACCS_keys_undecanes=data5[0].tolist()
ordered_MACCS_keys_undecanes=sorted(MACCS_keys_undecanes)
Morgan_circular_undecanes=data6[0].tolist()
ordered_Morgan_circular_undecanes=sorted(Morgan_circular_undecanes)


number=[]
for i in range (0,12561):
    number.append(i)
    
# Plotting
ax = plt.subplot(111)
plt.plot(number, ordered_adjacency_matrix_undecanes, label="Based on adjacency matrix ") 
plt.plot(number, ordered_atom_pair_undecanes, label="Based on atom pair")
plt.plot(number, ordered_topological_torsion_undecanes, label="Based on topological torsion") 
plt.plot(number, ordered_Avalon_undecanes, label="Based on Avalon") 
plt.plot(number, ordered_MACCS_keys_undecanes, label="Based on MACCS keys") 
plt.plot(number, ordered_Morgan_circular_undecanes, label="Based on Morgan circular") 
ax.legend(bbox_to_anchor=(1.1, 1.05)) 
plt.xlabel('Number of graphs to be compared')
plt.ylabel('Tanimoto coefficient')
plt.suptitle('For Undecanes')