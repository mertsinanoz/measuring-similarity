# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 14:27:20 2025

@author: sinan.oz
"""

import pandas as pd 
import matplotlib.pyplot as plt

data1=pd.read_excel('topological_indices_undecanes.xlsx')
data2=pd.read_excel('atom_pair_undecanes.xlsx')
data3=pd.read_excel('topological_torsion_undecanes.xlsx')
data4=pd.read_excel('Avalon_undecanes.xlsx')
data5=pd.read_excel('MACCS_keys_undecanes.xlsx')
data6=pd.read_excel('Morgan_circular_undecanes.xlsx')

# Relationships between Tanimoto coefficient values based on adjacency matrices
#and Tanimoto coefficient values based on molecular fingerprints for undecanes.
sorted_Topological_indices_undecanes=data1.sort_values(0)
sorted_atom_pair_undecanes=data2.sort_values(0)
sorted_topological_torsion_undecanes=data3.sort_values(0)
sorted_Avalon_undecanes=data4.sort_values(0)
sorted_MACCS_keys_undecanes=data5.sort_values(0)
sorted_Morgan_circular_undecanes=data6.sort_values(0)

# Plotting
ax = plt.subplot(111)
plt.plot(sorted_Topological_indices_undecanes,sorted_atom_pair_undecanes, linestyle='-', label="Topological indices vs atom pair")
plt.plot(sorted_Topological_indices_undecanes,sorted_topological_torsion_undecanes, linestyle='--', label="Topological indices vs topological torsion")
plt.plot(sorted_Topological_indices_undecanes,sorted_Avalon_undecanes, linestyle='-.', label="Topological indices vs Avalon")
plt.plot(sorted_Topological_indices_undecanes,sorted_MACCS_keys_undecanes, linestyle=':', label="Topological indices vs MACCS keys")
plt.plot(sorted_Topological_indices_undecanes,sorted_Morgan_circular_undecanes, linestyle='-', linewidth=2, label="Topological indices vs Morgan circular")

ax.legend(bbox_to_anchor=(1.1, 1.05)) 
plt.xlabel('Jaccard/Tanimoto index based on topological indices')
plt.ylabel('Jaccard/Tanimoto index based on molecular fingerprint')
plt.suptitle('For Undecanes')