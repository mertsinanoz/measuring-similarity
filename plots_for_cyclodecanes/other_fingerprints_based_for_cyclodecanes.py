# -*- coding: utf-8 -*-
"""
Created on Thu Jun  5 23:13:56 2025

@author: sinan.oz
"""

import pandas as pd 
import matplotlib.pyplot as plt

data2=pd.read_excel('atom_pair_cyclodecanes.xlsx')
data3=pd.read_excel('topological_torsion_cyclodecanes.xlsx')
data4=pd.read_excel('Avalon_cyclodecanes.xlsx')
data5=pd.read_excel('MACCS_keys_cyclodecanes.xlsx')
data6=pd.read_excel('Morgan_circular_cyclodecanes.xlsx')

#Relationships between Tanimoto coefficient values based on molecular fingerprints 
#for cyclodecanes
sorted_atom_pair_cyclodecanes=data2.sort_values(0)
sorted_topological_torsion_cyclodecanes=data3.sort_values(0)
sorted_Avalon_cyclodecanes=data4.sort_values(0)
sorted_MACCS_keys_cyclodecanes=data5.sort_values(0)
sorted_Morgan_circular_cyclodecanes=data6.sort_values(0)

# Plotting
ax = plt.subplot(111)
plt.plot(sorted_atom_pair_cyclodecanes,sorted_topological_torsion_cyclodecanes, label="Atom pair vs topological torsion")
plt.plot(sorted_atom_pair_cyclodecanes,sorted_Avalon_cyclodecanes, label="Atom pair vs Avalon")
plt.plot(sorted_atom_pair_cyclodecanes,sorted_MACCS_keys_cyclodecanes, label="Atom pair vs MACCS keys")
plt.plot(sorted_atom_pair_cyclodecanes,sorted_Morgan_circular_cyclodecanes, label="Atom pair vs Morgan circular")
plt.plot(sorted_topological_torsion_cyclodecanes,sorted_Avalon_cyclodecanes, label="Topological torsion vs Avalon")
plt.plot(sorted_topological_torsion_cyclodecanes,sorted_MACCS_keys_cyclodecanes, label="Topological torsion vs MACCS keys")
plt.plot(sorted_topological_torsion_cyclodecanes,sorted_Morgan_circular_cyclodecanes, label="Topological torsion vs Morgan circular")
plt.plot(sorted_Avalon_cyclodecanes,sorted_MACCS_keys_cyclodecanes, label="Avalon vs MACCS keys")
plt.plot(sorted_Avalon_cyclodecanes,sorted_Morgan_circular_cyclodecanes, label="Avalon vs Morgan circular")
plt.plot(sorted_MACCS_keys_cyclodecanes,sorted_Morgan_circular_cyclodecanes, label="MACCS keys vs Morgan circular")

ax.legend(bbox_to_anchor=(1.1, 1.05)) 
plt.xlabel('Jaccard/Tanimoto index')
plt.ylabel('Jaccard/Tanimoto index')
plt.suptitle('For Cyclodecanes')