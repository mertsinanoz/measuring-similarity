# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 14:35:46 2025

@author: sinan.oz
"""

import pandas as pd 
import matplotlib.pyplot as plt


data2=pd.read_excel('atom_pair_undecanes.xlsx')
data3=pd.read_excel('topological_torsion_undecanes.xlsx')
data4=pd.read_excel('Avalon_undecanes.xlsx')
data5=pd.read_excel('MACCS_keys_undecanes.xlsx')
data6=pd.read_excel('Morgan_circular_undecanes.xlsx')

#Relationships between Tanimoto coefficient values based on molecular fingerprints 
#for undecanes
sorted_atom_pair_undecanes=data2.sort_values(0)
sorted_topological_torsion_undecanes=data3.sort_values(0)
sorted_Avalon_undecanes=data4.sort_values(0)
sorted_MACCS_keys_undecanes=data5.sort_values(0)
sorted_Morgan_circular_undecanes=data6.sort_values(0)

# Plotting
ax = plt.subplot(111)
plt.plot(sorted_atom_pair_undecanes,sorted_topological_torsion_undecanes, label="Atom pair vs topological torsion")
plt.plot(sorted_atom_pair_undecanes,sorted_Avalon_undecanes, label="Atom pair vs Avalon")
plt.plot(sorted_atom_pair_undecanes,sorted_MACCS_keys_undecanes, label="Atom pair vs MACCS keys")
plt.plot(sorted_atom_pair_undecanes,sorted_Morgan_circular_undecanes, label="Atom pair vs Morgan circular")
plt.plot(sorted_topological_torsion_undecanes,sorted_Avalon_undecanes, label="Topological torsion vs Avalon")
plt.plot(sorted_topological_torsion_undecanes,sorted_MACCS_keys_undecanes, label="Topological torsion vs MACCS keys")
plt.plot(sorted_topological_torsion_undecanes,sorted_Morgan_circular_undecanes, label="Topological torsion vs Morgan circular")
plt.plot(sorted_Avalon_undecanes,sorted_MACCS_keys_undecanes, label="Avalon vs MACCS keys")
plt.plot(sorted_Avalon_undecanes,sorted_Morgan_circular_undecanes, label="Avalon vs Morgan circular")
plt.plot(sorted_MACCS_keys_undecanes,sorted_Morgan_circular_undecanes, label="MACCS keys vs Morgan circular")

ax.legend(bbox_to_anchor=(1.1, 1.05)) 
plt.xlabel('Jaccard/Tanimoto index')
plt.ylabel('Jaccard/Tanimoto index')
plt.suptitle('For Undecanes')