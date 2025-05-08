# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 15:33:35 2025

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


#Tanimoto coefficient values versus the number of tridecanes.
Topological_indices_tridecanes=data1[0].tolist()
ordered_Topological_indices_tridecanes=sorted(Topological_indices_tridecanes)
atom_pair_tridecanes=data2[0].tolist()
ordered_atom_pair_tridecanes=sorted(atom_pair_tridecanes)
topological_torsion_tridecanes=data3[0].tolist()
ordered_topological_torsion_tridecanes=sorted(topological_torsion_tridecanes)
Avalon_tridecanes=data4[0].tolist()
ordered_Avalon_tridecanes=sorted(Avalon_tridecanes)
MACCS_keys_tridecanes=data5[0].tolist()
ordered_MACCS_keys_tridecanes=sorted(MACCS_keys_tridecanes)
Morgan_circular_tridecanes=data6[0].tolist()
ordered_Morgan_circular_tridecanes=sorted(Morgan_circular_tridecanes)


number=[]
for i in range (0,321201):
    number.append(i)
 
# Plotting
ax = plt.subplot(111)
plt.plot(number, ordered_Topological_indices_tridecanes, label="Based on topological indices") 
plt.plot(number, ordered_atom_pair_tridecanes, label="Based on atom pair")
plt.plot(number, ordered_topological_torsion_tridecanes, label="Based on topological torsion") 
plt.plot(number, ordered_Avalon_tridecanes, label="Based on Avalon") 
plt.plot(number, ordered_MACCS_keys_tridecanes, label="Based on MACCS keys") 
plt.plot(number, ordered_Morgan_circular_tridecanes, label="Based on Morgan circular") 
ax.legend(bbox_to_anchor=(1.1, 1.05)) 
plt.xlabel('Number of graphs to be compared')
plt.ylabel('Jaccard/Tanimoto index')
plt.suptitle('For Tridecanes')
