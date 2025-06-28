# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 15:22:07 2025

@author: sinan.oz
"""

import pandas as pd 
import matplotlib.pyplot as plt

data1=pd.read_excel('topological_indices_dodecanes.xlsx')
data2=pd.read_excel('atom_pair_dodecanes.xlsx')
data3=pd.read_excel('topological_torsion_dodecanes.xlsx')
data4=pd.read_excel('Avalon_dodecanes.xlsx')
data5=pd.read_excel('MACCS_keys_dodecanes.xlsx')
data6=pd.read_excel('Morgan_circular_dodecanes.xlsx')

# Relationships between Tanimoto coefficient values based on adjacency matrices
#and Tanimoto coefficient values based on molecular fingerprints for dodecanes.
sorted_Topological_indices_dodecanes=data1.sort_values(0)
sorted_atom_pair_dodecanes=data2.sort_values(0)
sorted_topological_torsion_dodecanes=data3.sort_values(0)
sorted_Avalon_dodecanes=data4.sort_values(0)
sorted_MACCS_keys_dodecanes=data5.sort_values(0)
sorted_Morgan_circular_dodecanes=data6.sort_values(0)

# Plotting
ax = plt.subplot(111)
plt.plot(sorted_Topological_indices_dodecanes,sorted_atom_pair_dodecanes, linestyle='-', label="Topological indices vs atom pair ")
plt.plot(sorted_Topological_indices_dodecanes,sorted_topological_torsion_dodecanes, linestyle='--', label="Topological indices vs topological torsion")
plt.plot(sorted_Topological_indices_dodecanes,sorted_Avalon_dodecanes, linestyle='-.', label="Topological indices vs Avalon")
plt.plot(sorted_Topological_indices_dodecanes,sorted_MACCS_keys_dodecanes, linestyle=':', label="Topological indices vs MACCS keys")
plt.plot(sorted_Topological_indices_dodecanes,sorted_Morgan_circular_dodecanes, linestyle='-', linewidth=2, label="Topological indices vs Morgan circular")

ax.legend(bbox_to_anchor=(1.1, 1.05)) 
plt.xlabel('Jaccard/Tanimoto index based on topological indices')
plt.ylabel('Jaccard/Tanimoto index based on molecular fingerprint')
plt.suptitle('For Dodecanes')