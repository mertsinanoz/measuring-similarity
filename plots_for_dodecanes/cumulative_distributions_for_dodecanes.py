# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 15:23:13 2025

@author: sinan.oz
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np 

df1=pd.read_excel('topological_indices_dodecanes.xlsx')
df2=pd.read_excel('atom_pair_dodecanes.xlsx')
df3=pd.read_excel('topological_torsion_dodecanes.xlsx')
df4=pd.read_excel('Avalon_dodecanes.xlsx')
df5=pd.read_excel('MACCS_keys_dodecanes.xlsx')
df6=pd.read_excel('Morgan_circular_dodecanes.xlsx')


#Cumulative distributions of the Tanimoto coefficients based on 
#the molecular fingerprints and adjacency matrices of dodecanes.
data_1 = df1.iloc[:, 0].dropna()
data_2 = df2.iloc[:, 0].dropna()
data_3 = df3.iloc[:, 0].dropna()
data_4 = df4.iloc[:, 0].dropna()
data_5 = df5.iloc[:, 0].dropna()
data_6 = df6.iloc[:, 0].dropna()

sorted_data1 = np.sort(data_1)
sorted_data2 = np.sort(data_2)
sorted_data3 = np.sort(data_3)
sorted_data4 = np.sort(data_4)
sorted_data5 = np.sort(data_5)
sorted_data6 = np.sort(data_6)

cdf1= np.cumsum(data_1) / np.sum(data_1)
cdf2= np.cumsum(data_2) / np.sum(data_2)
cdf3= np.cumsum(data_3) / np.sum(data_3)
cdf4= np.cumsum(data_4) / np.sum(data_4)
cdf5= np.cumsum(data_5) / np.sum(data_5)
cdf6= np.cumsum(data_6) / np.sum(data_6)

# Plotting
ax = plt.subplot(111)
plt.plot(sorted_data1,cdf1,label="Based on topological indices")
plt.plot(sorted_data2,cdf2,label="Based on atom pair")
plt.plot(sorted_data3,cdf3,label="Based on topological torsion")
plt.plot(sorted_data4,cdf4,label="Based on Avalon")
plt.plot(sorted_data5,cdf5,label="Based on MACCS keys")
plt.plot(sorted_data6,cdf6,label="Based on Morgan circular")
ax.legend(bbox_to_anchor=(1.1, 1.05))
plt.xlabel('Jaccard/Tanimoto index')
plt.ylabel('%')
plt.suptitle('For Dodecanes')