# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 14:27:20 2025

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

# Relationships between Tanimoto coefficient values based on adjacency matrices
#and Tanimoto coefficient values based on molecular fingerprints for undecanes.
sorted_adjacency_matrix_undecanes=data1.sort_values(0)
sorted_atom_pair_undecanes=data2.sort_values(0)
sorted_topological_torsion_undecanes=data3.sort_values(0)
sorted_Avalon_undecanes=data4.sort_values(0)
sorted_MACCS_keys_undecanes=data5.sort_values(0)
sorted_Morgan_circular_undecanes=data6.sort_values(0)

# Plotting
ax = plt.subplot(111)
plt.plot(sorted_adjacency_matrix_undecanes,sorted_atom_pair_undecanes, label="Adjacency matrix vs atom pair ")
plt.plot(sorted_adjacency_matrix_undecanes,sorted_topological_torsion_undecanes, label="Adjacency matrix vs topological torsion")
plt.plot(sorted_adjacency_matrix_undecanes,sorted_Avalon_undecanes, label="Adjacency matrix vs Avalon")
plt.plot(sorted_adjacency_matrix_undecanes,sorted_MACCS_keys_undecanes, label="Adjacency matrix vs MACCS keys")
plt.plot(sorted_adjacency_matrix_undecanes,sorted_Morgan_circular_undecanes, label="Adjacency matrix vs Morgan circular")

ax.legend(bbox_to_anchor=(1.1, 1.05)) 
plt.xlabel('Tanimoto coefficient based on adjacency matrix')
plt.ylabel('Tanimoto coefficient based on molecular fingerprint')
plt.suptitle('For Undecanes')