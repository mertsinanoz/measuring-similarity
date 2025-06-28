# -*- coding: utf-8 -*-
"""
Created on Thu Jun 12 23:13:42 2025

@author: sinan.oz
"""

import cirpy as cp
import numpy as np
from pysmiles import read_smiles
import networkx as nx
import pickle
import pandas as pd
import math
from typing import Set
import matplotlib.pyplot as plt

with open("./isomer11.pkl", "rb") as file11:
    isomer11= pickle.load(file11)

#Sorting the graphs in the isomer11 file according to the SMILES order in the isomers_of_undecanes file
df=pd.read_excel('isomers_of_undecanes.xlsx')
graphs_of_undecanes=[None] * 159
for i in range(len(df)):
   smi=cp.resolve(df.Isomer[i], 'smiles')
   G = read_smiles(smi)
   for j in range(len(isomer11)):
       is_isomorphic=nx.is_isomorphic(G, isomer11[j])
       if is_isomorphic:
         graphs_of_undecanes[i]=G

def tanimoto_coefficient(set_1: Set[int], set_2: Set[int]) -> float:
   intersect= len(set_1.intersection(set_2))
   union = len(set_1) + len(set_2) - intersect
   return intersect / union

def Degree_based_indices(G):
   edgeTable1 = G.edges
   Q = list(edgeTable1)
   sum1= 0
   sum2=0
   sum3=0
   sum4=0
   sum5= 0
   sum6=0
   sum7=0
   sum8=0
   sum9=0
   for i in range(len(Q)):
      j= 0
      #Sigma Index
      k1 = (G.degree(Q[i][j])-G.degree(Q[i][j+1]))**2
      sum1 += k1

      #The Albertson index
      k2 = abs(G.degree(Q[i][j])-G.degree(Q[i][j+1]))
      sum2 += k2
        
      #The GA index
      k3 = (2*math.sqrt(G.degree(Q[i][j])*G.degree(Q[i][j+1]))/
            (G.degree(Q[i][j])+G.degree(Q[i][j+1])))
      sum3 += k3
        
      #The Sum connectivity
      k4 = 1/math.sqrt(G.degree(Q[i][j])+G.degree(Q[i][j+1]))
      sum4 += k4
        
      #The ISI index
      k5 = ((G.degree(Q[i][j])*G.degree(Q[i][j+1]))/
            (G.degree(Q[i][j])+G.degree(Q[i][j+1])))
      sum5 += k5
    
      #ABC
      k6 = math.sqrt((G.degree(Q[i][j])+G.degree(Q[i][j+1])-2)/
          (G.degree(Q[i][j])*G.degree(Q[i][j+1])))
      sum6 += k6
      
      #AZI index
      k7 = ((G.degree(Q[i][j])*G.degree(Q[i][j+1]))/
            (G.degree(Q[i][j])+G.degree(Q[i][j+1])-2))**3
      sum7 += k7
      
      #The SDD index
      k8 = ((G.degree(Q[i][j])/G.degree(Q[i][j+1]))+
      (G.degree(Q[i][j+1])/G.degree(Q[i][j])))
      sum8 += k8
      
      #The Randic index
      k9= 1/math.sqrt(G.degree(Q[i][j])*G.degree(Q[i][j+1]))
      sum9 += k9     
   return sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,sum9

def Wiener_and_Balaban_indices(G):
   m= G.number_of_edges()
   n= G.number_of_nodes()
   shortest_paths = dict(nx.all_pairs_shortest_path_length(G))
   node_list = list(G.nodes())
   distance_matrix = np.zeros((len(node_list), len(node_list)))
   for i, node1 in enumerate(node_list):
      for j, node2 in enumerate(node_list):
         distance_matrix[i, j] = shortest_paths[node1].get(node2, np.inf)
   sum_of_entries=np.sum(distance_matrix)      
   matrix=np.sum(distance_matrix, axis=1)
   row_sums=[]
   for i in range(len(matrix)):
      row_sums.append(matrix[i])
   J = 0
   edgeTable1 = G.edges
   Q1 = list(edgeTable1) 
   for index1 in range(len(Q1)):
      J +=1/math.sqrt(row_sums[Q1[index1][0]-1]*row_sums[Q1[index1][1]-1])   
   return sum_of_entries/2, m*J/(m-n+2)

def Mostar_Szeged_index(G):
   vertices=list(G.nodes)
   edge_table = G.edges
   Q = list(edge_table)
   totalMostar=0
   totalSzeged=0
   for i1 in range(len(Q)):
       j1=0
       nu1=0
       nv1=0
       for p in range(len(vertices)):
           if len(nx.shortest_path(G, Q[i1][j1], p))<len(nx.shortest_path(G, Q[i1][j1+1], p)):
               nu1=nu1+1
           if len(nx.shortest_path(G, Q[i1][j1], p))>len(nx.shortest_path(G, Q[i1][j1+1], p)):
               nv1=nv1+1
       totalMostar += abs(nu1-nv1)
       totalSzeged += nu1*nv1
   return totalMostar,totalSzeged

def variance_of_graph(G):
   n=G.number_of_nodes()
   sum1= 0
   sum2=0
   for i in range(n):
      k1=(G.degree[i]**2)
      sum1 += k1
      k2=G.degree[i]
      sum2 += k2
   variance=1/n*sum1-(1/n*sum2)**2
   return variance

def total_irregularity(G):
   total_irregularity= 0
   for i in range(0,G.number_of_nodes()-1):
      for j in range(i+1,G.number_of_nodes()):
         k1=abs(G.degree[i]-G.degree[j])
         total_irregularity += k1
   return total_irregularity

def energy_and_spectral_radius_of_graph(G):
   A = nx.adjacency_matrix(G)
   eigenvalues = np.linalg.eigvals(A.toarray())
   spectral_radius = max(abs(eigenvalues))
   summation=0
   for i in range(len(eigenvalues)):
      summation+=round(abs(eigenvalues[i]),5)
   return summation,spectral_radius

def centrality_measures(G):
    degree_centrality = nx.degree_centrality(G)
    degree_values_list = list(degree_centrality.values())
    closeness_centrality = nx.closeness_centrality(G)
    closeness_values_list = list(closeness_centrality.values())
    betweenness_centrality = nx.betweenness_centrality(G, normalized = True, endpoints = False)
    betweenness_values_list = list(betweenness_centrality .values())
    degree_index=sum(degree_values_list)
    closeness_index=sum(closeness_values_list)
    betweenness_index=sum(betweenness_values_list)
    return degree_index,closeness_index,betweenness_index

#Tanimoto coefficient calculations based on topological indices for undecanes 
#and exporting the obtained data as pkl and xlsx files.
topological_indices_undecanes=[]
def similarity(vector1, vector2, threshold=0.5):
    difference = np.abs(np.array(vector1) - np.array(vector2))
    intersection = 0
    union = 0
    for d in difference:
        if d < threshold:
            intersection += 1
        else:
            union += 1
    union += intersection
    return intersection / union

thresholds = np.linspace(0, 1, 10)   # 10 different threshold values between 0 and 1
num_graphs = len(graphs_of_undecanes)
all_similarity_results = []
for i in range(num_graphs-1):
   for j in range(i+1, num_graphs):
      vector1=[]
      vector2=[]
      degree_based_indices1=Degree_based_indices(graphs_of_undecanes[i])
      centrality_measures1=centrality_measures(graphs_of_undecanes[i])
      energy_and_spectral_radius1=energy_and_spectral_radius_of_graph(graphs_of_undecanes[i])
      vector1.append(round(degree_based_indices1[0],5))
      vector1.append(round(degree_based_indices1[1],5))
      vector1.append(round(degree_based_indices1[2],5))
      vector1.append(round(degree_based_indices1[3],5))
      vector1.append(round(degree_based_indices1[4],5))
      vector1.append(round(degree_based_indices1[5],5))
      vector1.append(round(degree_based_indices1[6],5))
      vector1.append(round(degree_based_indices1[7],5))
      vector1.append(round(degree_based_indices1[8],5))
      vector1.append(round(Mostar_Szeged_index(graphs_of_undecanes[i])[0],5))
      vector1.append(round(Mostar_Szeged_index(graphs_of_undecanes[i])[1],5))
      vector1.append(round(variance_of_graph(graphs_of_undecanes[i]),5))
      vector1.append(round(total_irregularity(graphs_of_undecanes[i]),5))
      vector1.append(round(Wiener_and_Balaban_indices(graphs_of_undecanes[i])[0],5))
      vector1.append(round(Wiener_and_Balaban_indices(graphs_of_undecanes[i])[1],5))
      vector1.append(round(energy_and_spectral_radius1[0],5))
      vector1.append(round(energy_and_spectral_radius1[1],5))
      vector1.append(round(centrality_measures1[0],5))
      vector1.append(round(centrality_measures1[1],5))
      vector1.append(round(centrality_measures1[2],5))

    
      degree_based_indices2=Degree_based_indices(graphs_of_undecanes[j])
      centrality_measures2=centrality_measures(graphs_of_undecanes[j])
      energy_and_spectral_radius2=energy_and_spectral_radius_of_graph(graphs_of_undecanes[j])
      vector2.append(round(degree_based_indices2[0],5))
      vector2.append(round(degree_based_indices2[1],5))
      vector2.append(round(degree_based_indices2[2],5))
      vector2.append(round(degree_based_indices2[3],5))
      vector2.append(round(degree_based_indices2[4],5))
      vector2.append(round(degree_based_indices2[5],5))
      vector2.append(round(degree_based_indices2[6],5))
      vector2.append(round(degree_based_indices2[7],5))
      vector2.append(round(degree_based_indices2[8],5))
      vector2.append(round(Mostar_Szeged_index(graphs_of_undecanes[j])[0],5))
      vector2.append(round(Mostar_Szeged_index(graphs_of_undecanes[j])[1],5))
      vector2.append(round(variance_of_graph(graphs_of_undecanes[j]),5))
      vector2.append(round(total_irregularity(graphs_of_undecanes[j]),5))
      vector2.append(round(Wiener_and_Balaban_indices(graphs_of_undecanes[j])[0],5))
      vector2.append(round(Wiener_and_Balaban_indices(graphs_of_undecanes[j])[1],5))
      vector2.append(round(energy_and_spectral_radius2[0],5))
      vector2.append(round(energy_and_spectral_radius2[1],5))
      vector2.append(round(centrality_measures2[0],5))
      vector2.append(round(centrality_measures2[1],5))
      vector2.append(round(centrality_measures2[2],5))
      
      # Similarity calculation for different thresholds
      similarity_for_current_pair = []
      for threshold in thresholds:
         similarity_for_current_pair.append(similarity(vector1, vector2, threshold))
      all_similarity_results.append(similarity_for_current_pair)

for i, similarity_values in enumerate(all_similarity_results[::106]):  
    plt.plot(thresholds, similarity_values, label=f'Graph Pair {i+1}')
plt.xlabel("Threshold"), plt.ylabel("Similarity")
plt.title("Sensitivity Analysis on the selection of 0.5")
plt.show()