# -*- coding: utf-8 -*-
"""
Created on Thu Jun  5 23:37:54 2025

@author: sinan.oz
"""

import pandas as pd 
import matplotlib.pyplot as plt

data1=pd.read_excel('topological_indices_connected_graphs_with_7_vertices.xlsx')

#Tanimoto coefficient values versus the number of undecanes.
Topological_indices_undecanes=data1[0].tolist()
ordered_Topological_indices_undecanes=sorted(Topological_indices_undecanes)

number=[]
for i in range (0,363378):
    number.append(i)

# Plotting
ax = plt.subplot(111)
plt.plot(number, ordered_Topological_indices_undecanes, label="Based on topological indices ") 
ax.legend(bbox_to_anchor=(1.1, 1.05)) 
plt.xlabel('Number of graphs to be compared')
plt.ylabel('Jaccard/Tanimoto index')
plt.suptitle('For Connected Graphs with 7 Vertices')