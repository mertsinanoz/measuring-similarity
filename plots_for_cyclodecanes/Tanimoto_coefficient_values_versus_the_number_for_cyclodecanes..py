# -*- coding: utf-8 -*-
"""
Created on Thu Jun  5 23:15:52 2025

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

#Tanimoto coefficient values versus the number of cyclodecanes.
Topological_indices_cyclodecanes=data1[0].tolist()
ordered_Topological_indices_cyclodecanes=sorted(Topological_indices_cyclodecanes)
atom_pair_cyclodecanes=data2[0].tolist()
ordered_atom_pair_cyclodecanes=sorted(atom_pair_cyclodecanes)
topological_torsion_cyclodecanes=data3[0].tolist()
ordered_topological_torsion_cyclodecanes=sorted(topological_torsion_cyclodecanes)
Avalon_cyclodecanes=data4[0].tolist()
ordered_Avalon_cyclodecanes=sorted(Avalon_cyclodecanes)
MACCS_keys_cyclodecanes=data5[0].tolist()
ordered_MACCS_keys_cyclodecanes=sorted(MACCS_keys_cyclodecanes)
Morgan_circular_cyclodecanes=data6[0].tolist()
ordered_Morgan_circular_cyclodecanes=sorted(Morgan_circular_cyclodecanes)

number=[]
for i in range (0,112575):
    number.append(i)
    
# Plotting
ax = plt.subplot(111)
plt.plot(number, ordered_Topological_indices_cyclodecanes, linestyle='-', label="Based on topological indices") 
plt.plot(number, ordered_atom_pair_cyclodecanes, linestyle='--', label="Based on atom pair")
plt.plot(number, ordered_topological_torsion_cyclodecanes, linestyle='-.', label="Based on topological torsion") 
plt.plot(number, ordered_Avalon_cyclodecanes, linestyle=':', label="Based on Avalon") 
plt.plot(number, ordered_MACCS_keys_cyclodecanes, linestyle='-', linewidth=0.5, label="Based on MACCS keys") 
plt.plot(number, ordered_Morgan_circular_cyclodecanes, linestyle='-', linewidth=2, label="Based on Morgan circular") 
ax.legend(bbox_to_anchor=(1.1, 1.05))
plt.xlabel('Number of graphs to be compared')
plt.ylabel('Jaccard/Tanimoto index')
plt.suptitle('For Cyclodecanes')