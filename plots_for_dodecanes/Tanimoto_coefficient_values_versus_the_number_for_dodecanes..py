# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 15:16:57 2025

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

#Tanimoto coefficient values versus the number of dodecanes.
Topological_indices_dodecanes=data1[0].tolist()
ordered_Topological_indices_dodecanes=sorted(Topological_indices_dodecanes)
atom_pair_dodecanes=data2[0].tolist()
ordered_atom_pair_dodecanes=sorted(atom_pair_dodecanes)
topological_torsion_dodecanes=data3[0].tolist()
ordered_topological_torsion_dodecanes=sorted(topological_torsion_dodecanes)
Avalon_dodecanes=data4[0].tolist()
ordered_Avalon_dodecanes=sorted(Avalon_dodecanes)
MACCS_keys_dodecanes=data5[0].tolist()
ordered_MACCS_keys_dodecanes=sorted(MACCS_keys_dodecanes)
Morgan_circular_dodecanes=data6[0].tolist()
ordered_Morgan_circular_dodecanes=sorted(Morgan_circular_dodecanes)

number=[]
for i in range (0,62835):
    number.append(i)
    
# Plotting
ax = plt.subplot(111)
plt.plot(number, ordered_Topological_indices_dodecanes, linestyle='-', label="Based on topological indices") 
plt.plot(number, ordered_atom_pair_dodecanes, linestyle='--', label="Based on atom pair")
plt.plot(number, ordered_topological_torsion_dodecanes, linestyle='-.', label="Based on topological torsion") 
plt.plot(number, ordered_Avalon_dodecanes, linestyle=':', label="Based on Avalon") 
plt.plot(number, ordered_MACCS_keys_dodecanes, linestyle='-', linewidth=0.5, label="Based on MACCS keys") 
plt.plot(number, ordered_Morgan_circular_dodecanes, linestyle='-', linewidth=2, label="Based on Morgan circular") 
ax.legend(bbox_to_anchor=(1.1, 1.05)) 
plt.xlabel('Number of graphs to be compared')
plt.ylabel('Jaccard/Tanimoto index')
plt.suptitle('For Dodecanes')