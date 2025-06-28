# -*- coding: utf-8 -*-
"""
Created on Mon Jun  9 21:39:37 2025

@author: sinan.oz
"""

import pandas as pd 
import matplotlib.pyplot as plt

data2=pd.read_excel('atom_pair_cycloundecanes.xlsx')
data3=pd.read_excel('topological_torsion_cycloundecanes.xlsx')
data4=pd.read_excel('Avalon_cycloundecanes.xlsx')
data5=pd.read_excel('MACCS_keys_cycloundecanes.xlsx')
data6=pd.read_excel('Morgan_circular_cycloundecanes.xlsx')

#Relationships between Tanimoto coefficient values based on molecular fingerprints 
#for cycloundecanes
sorted_atom_pair_cycloundecanes=data2.sort_values(0)
sorted_topological_torsion_cycloundecanes=data3.sort_values(0)
sorted_Avalon_cycloundecanes=data4.sort_values(0)
sorted_MACCS_keys_cycloundecanes=data5.sort_values(0)
sorted_Morgan_circular_cycloundecanes=data6.sort_values(0)

# Plotting
ax = plt.subplot(111)
plt.plot(sorted_atom_pair_cycloundecanes,sorted_topological_torsion_cycloundecanes, label="Atom pair vs topological torsion")
plt.plot(sorted_atom_pair_cycloundecanes,sorted_Avalon_cycloundecanes, label="Atom pair vs Avalon")
plt.plot(sorted_atom_pair_cycloundecanes,sorted_MACCS_keys_cycloundecanes, label="Atom pair vs MACCS keys")
plt.plot(sorted_atom_pair_cycloundecanes,sorted_Morgan_circular_cycloundecanes, label="Atom pair vs Morgan circular")
plt.plot(sorted_topological_torsion_cycloundecanes,sorted_Avalon_cycloundecanes, label="Topological torsion vs Avalon")
plt.plot(sorted_topological_torsion_cycloundecanes,sorted_MACCS_keys_cycloundecanes, label="Topological torsion vs MACCS keys")
plt.plot(sorted_topological_torsion_cycloundecanes,sorted_Morgan_circular_cycloundecanes, label="Topological torsion vs Morgan circular")
plt.plot(sorted_Avalon_cycloundecanes,sorted_MACCS_keys_cycloundecanes, label="Avalon vs MACCS keys")
plt.plot(sorted_Avalon_cycloundecanes,sorted_Morgan_circular_cycloundecanes, label="Avalon vs Morgan circular")
plt.plot(sorted_MACCS_keys_cycloundecanes,sorted_Morgan_circular_cycloundecanes, label="MACCS keys vs Morgan circular")

ax.legend(bbox_to_anchor=(1.1, 1.05)) 
plt.xlabel('Jaccard/Tanimoto index')
plt.ylabel('Jaccard/Tanimoto index')
plt.suptitle('For Cycloundecanes')