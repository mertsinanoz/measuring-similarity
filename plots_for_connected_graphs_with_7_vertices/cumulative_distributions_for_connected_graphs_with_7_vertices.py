# -*- coding: utf-8 -*-
"""
Created on Thu Jun  5 23:34:23 2025

@author: sinan.oz
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np 

df1=pd.read_excel('topological_indices_connected_graphs_with_7_vertices.xlsx')

#Cumulative distributions of the Tanimoto coefficients based on 
#the molecular fingerprints and adjacency matrices of cyclodecanes.
data_1 = df1.iloc[:, 0].dropna()
sorted_data1 = np.sort(data_1)
cdf1= np.cumsum(data_1) / np.sum(data_1)

# Plotting
ax = plt.subplot(111)
plt.plot(sorted_data1,cdf1,label="Based on topological indices")
ax.legend(bbox_to_anchor=(1.1, 1.05))
plt.xlabel('Jaccard/Tanimoto index')
plt.ylabel('%')
plt.suptitle('For Connected Graphs with 7 Vertices')