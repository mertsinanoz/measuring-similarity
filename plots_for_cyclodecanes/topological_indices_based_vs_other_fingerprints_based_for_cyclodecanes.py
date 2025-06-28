# -*- coding: utf-8 -*-
"""
Created on Thu Jun  5 23:29:07 2025

@author: sinan.oz
"""

import pandas as pd 
import matplotlib.pyplot as plt

data1=pd.read_excel('topological_indices_cyclodecanes.xlsx')
data2=pd.read_excel('atom_pair_cyclodecanes.xlsx')
data3=pd.read_excel('topological_torsion_cyclodecanes.xlsx')
data4=pd.read_excel('Avalon_cyclodecanes.xlsx')
data5=pd.read_excel('MACCS_keys_cyclodecanes.xlsx')
data6=pd.read_excel('Morgan_circular_cyclodecanes.xlsx')

# Relationships between Tanimoto coefficient values based on adjacency matrices
#and Tanimoto coefficient values based on molecular fingerprints for cyclodecanes.
sorted_Topological_indices_cyclodecanes=data1.sort_values(0)
sorted_atom_pair_cyclodecanes=data2.sort_values(0)
sorted_topological_torsion_cyclodecanes=data3.sort_values(0)
sorted_Avalon_cyclodecanes=data4.sort_values(0)
sorted_MACCS_keys_cyclodecanes=data5.sort_values(0)
sorted_Morgan_circular_cyclodecanes=data6.sort_values(0)

# Plotting
ax = plt.subplot(111)
plt.plot(sorted_Topological_indices_cyclodecanes,sorted_atom_pair_cyclodecanes, linestyle='-', label="Topological indices vs atom pair ")
plt.plot(sorted_Topological_indices_cyclodecanes,sorted_topological_torsion_cyclodecanes, linestyle='--', label="Topological indices vs topological torsion")
plt.plot(sorted_Topological_indices_cyclodecanes,sorted_Avalon_cyclodecanes, linestyle='-.', label="Topological indices vs Avalon")
plt.plot(sorted_Topological_indices_cyclodecanes,sorted_MACCS_keys_cyclodecanes, linestyle=':', label="Topological indices vs MACCS keys")
plt.plot(sorted_Topological_indices_cyclodecanes,sorted_Morgan_circular_cyclodecanes, linestyle='-', linewidth=2, label="Topological indices vs Morgan circular")

ax.legend(bbox_to_anchor=(1.1, 1.05)) 
plt.xlabel('Jaccard/Tanimoto index based on topological indices')
plt.ylabel('Jaccard/Tanimoto index based on molecular fingerprint')
plt.suptitle('For Cyclodecanes')