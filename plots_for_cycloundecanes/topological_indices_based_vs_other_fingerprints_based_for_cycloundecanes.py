# -*- coding: utf-8 -*-
"""
Created on Mon Jun  9 21:46:11 2025

@author: sinan.oz
"""

import pandas as pd 
import matplotlib.pyplot as plt

data1=pd.read_excel('topological_indices_cycloundecanes.xlsx')
data2=pd.read_excel('atom_pair_cycloundecanes.xlsx')
data3=pd.read_excel('topological_torsion_cycloundecanes.xlsx')
data4=pd.read_excel('Avalon_cycloundecanes.xlsx')
data5=pd.read_excel('MACCS_keys_cycloundecanes.xlsx')
data6=pd.read_excel('Morgan_circular_cycloundecanes.xlsx')

# Relationships between Tanimoto coefficient values based on adjacency matrices
#and Tanimoto coefficient values based on molecular fingerprints for cycloundecanes.
sorted_Topological_indices_cycloundecanes=data1.sort_values(0)
sorted_atom_pair_cycloundecanes=data2.sort_values(0)
sorted_topological_torsion_cycloundecanes=data3.sort_values(0)
sorted_Avalon_cycloundecanes=data4.sort_values(0)
sorted_MACCS_keys_cycloundecanes=data5.sort_values(0)
sorted_Morgan_circular_cycloundecanes=data6.sort_values(0)

# Plotting
ax = plt.subplot(111)
plt.plot(sorted_Topological_indices_cycloundecanes,sorted_atom_pair_cycloundecanes, linestyle='-', label="Topological indices vs atom pair ")
plt.plot(sorted_Topological_indices_cycloundecanes,sorted_topological_torsion_cycloundecanes, linestyle='--', label="Topological indices vs topological torsion")
plt.plot(sorted_Topological_indices_cycloundecanes,sorted_Avalon_cycloundecanes, linestyle='-.', label="Topological indices vs Avalon")
plt.plot(sorted_Topological_indices_cycloundecanes,sorted_MACCS_keys_cycloundecanes, linestyle=':', label="Topological indices vs MACCS keys")
plt.plot(sorted_Topological_indices_cycloundecanes,sorted_Morgan_circular_cycloundecanes, linestyle='-', linewidth=2, label="Topological indices vs Morgan circular")

ax.legend(bbox_to_anchor=(1.1, 1.05)) 
plt.xlabel('Jaccard/Tanimoto index based on topological indices')
plt.ylabel('Jaccard/Tanimoto index based on molecular fingerprint')
plt.suptitle('For Cycloundecanes')