# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 16:49:32 2025

@author: sinan.oz
"""

import pandas as pd 
import matplotlib.pyplot as plt

data1=pd.read_excel('topological_indices_tridecanes.xlsx')
data2=pd.read_excel('atom_pair_tridecanes.xlsx')
data3=pd.read_excel('topological_torsion_tridecanes.xlsx')
data4=pd.read_excel('Avalon_tridecanes.xlsx')
data5=pd.read_excel('MACCS_keys_tridecanes.xlsx')
data6=pd.read_excel('Morgan_circular_tridecanes.xlsx')

# Relationships between Tanimoto coefficient values based on adjacency matrices
#and Tanimoto coefficient values based on molecular fingerprints for tridecanes.
sorted_Topological_indices_tridecanes=data1.sort_values(0)
sorted_atom_pair_tridecanes=data2.sort_values(0)
sorted_topological_torsion_tridecanes=data3.sort_values(0)
sorted_Avalon_tridecanes=data4.sort_values(0)
sorted_MACCS_keys_tridecanes=data5.sort_values(0)
sorted_Morgan_circular_tridecanes=data6.sort_values(0)

# Plotting
ax = plt.subplot(111)
plt.plot(sorted_Topological_indices_tridecanes,sorted_atom_pair_tridecanes, label="Topological indices vs atom pair ")
plt.plot(sorted_Topological_indices_tridecanes,sorted_topological_torsion_tridecanes, label="Topological indices vs topological torsion")
plt.plot(sorted_Topological_indices_tridecanes,sorted_Avalon_tridecanes, label="Topological indices vs Avalon")
plt.plot(sorted_Topological_indices_tridecanes,sorted_MACCS_keys_tridecanes, label="Topological indices vs MACCS keys")
plt.plot(sorted_Topological_indices_tridecanes,sorted_Morgan_circular_tridecanes, label="Topological indices vs Morgan circular")

ax.legend(bbox_to_anchor=(1.1, 1.05)) 
plt.xlabel('Jaccard/Tanimoto index based on topological indices')
plt.ylabel('Jaccard/Tanimoto index based on molecular fingerprint')
plt.suptitle('For Tridecanes')