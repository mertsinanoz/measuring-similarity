# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 10:54:00 2025

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


#Tanimoto coefficient values versus the number of undecanes.
Topological_indices_undecanes=data1[0].tolist()
ordered_Topological_indices_undecanes=sorted(Topological_indices_undecanes)
atom_pair_undecanes=data2[0].tolist()
ordered_atom_pair_undecanes=sorted(atom_pair_undecanes)
topological_torsion_undecanes=data3[0].tolist()
ordered_topological_torsion_undecanes=sorted(topological_torsion_undecanes)
Avalon_undecanes=data4[0].tolist()
ordered_Avalon_undecanes=sorted(Avalon_undecanes)
MACCS_keys_undecanes=data5[0].tolist()
ordered_MACCS_keys_undecanes=sorted(MACCS_keys_undecanes)
Morgan_circular_undecanes=data6[0].tolist()
ordered_Morgan_circular_undecanes=sorted(Morgan_circular_undecanes)


number=[]
for i in range (0,12561):
    number.append(i)
    
# Plotting
ax = plt.subplot(111)
plt.plot(number, ordered_Topological_indices_undecanes, label="Based on topological indices ") 
plt.plot(number, ordered_atom_pair_undecanes, label="Based on atom pair")
plt.plot(number, ordered_topological_torsion_undecanes, label="Based on topological torsion") 
plt.plot(number, ordered_Avalon_undecanes, label="Based on Avalon") 
plt.plot(number, ordered_MACCS_keys_undecanes, label="Based on MACCS keys") 
plt.plot(number, ordered_Morgan_circular_undecanes, label="Based on Morgan circular") 
ax.legend(bbox_to_anchor=(1.1, 1.05)) 
plt.xlabel('Number of graphs to be compared')
plt.ylabel('Jaccard/Tanimoto index')
plt.suptitle('For Undecanes')