# -*- coding: utf-8 -*-
"""
Created on Mon Jun  9 21:42:35 2025

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

#Tanimoto coefficient values versus the number of cycloundecanes.
Topological_indices_cycloundecanes=data1[0].tolist()
ordered_Topological_indices_cycloundecanes=sorted(Topological_indices_cycloundecanes)
atom_pair_cycloundecanes=data2[0].tolist()
ordered_atom_pair_cycloundecanes=sorted(atom_pair_cycloundecanes)
topological_torsion_cycloundecanes=data3[0].tolist()
ordered_topological_torsion_cycloundecanes=sorted(topological_torsion_cycloundecanes)
Avalon_cycloundecanes=data4[0].tolist()
ordered_Avalon_cycloundecanes=sorted(Avalon_cycloundecanes)
MACCS_keys_cycloundecanes=data5[0].tolist()
ordered_MACCS_keys_cycloundecanes=sorted(MACCS_keys_cycloundecanes)
Morgan_circular_cycloundecanes=data6[0].tolist()
ordered_Morgan_circular_cycloundecanes=sorted(Morgan_circular_cycloundecanes)

number=[]
for i in range (0,757065):
    number.append(i)
    
# Plotting
ax = plt.subplot(111)
plt.plot(number, ordered_Topological_indices_cycloundecanes, linestyle='-', label="Based on topological indices") 
plt.plot(number, ordered_atom_pair_cycloundecanes, linestyle='--', label="Based on atom pair")
plt.plot(number, ordered_topological_torsion_cycloundecanes, linestyle='-.', label="Based on topological torsion") 
plt.plot(number, ordered_Avalon_cycloundecanes, linestyle=':', label="Based on Avalon") 
plt.plot(number, ordered_MACCS_keys_cycloundecanes, linestyle='-', linewidth=0.5, label="Based on MACCS keys") 
plt.plot(number, ordered_Morgan_circular_cycloundecanes, linestyle='-', linewidth=2, label="Based on Morgan circular") 
ax.legend(bbox_to_anchor=(1.1, 1.05))
plt.xlabel('Number of graphs to be compared')
plt.ylabel('Jaccard/Tanimoto index')
plt.suptitle('For Cycloundecanes')