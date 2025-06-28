# -*- coding: utf-8 -*-
"""
Created on Thu Jun  5 18:21:13 2025

@author: sinan.oz
"""

import numpy as np
import networkx as nx
import pickle
import pandas as pd
import math
from typing import Set
       
connected_graphsS7=[]   
with open('./data_acquisition_connected_graphs_with_7_vertices/connected_graphs_7_vertices.g6','r') as f:
     for line in f:
       connected_graphsS7.append(line.strip())
       
graphs_with_7_vertices=[]     
for i in range(len(connected_graphsS7)):
   str1=connected_graphsS7[i]
   g1=nx.from_graph6_bytes(bytes(str1, 'ascii', errors='strict'))
   graphs_with_7_vertices.append(g1) 

       
def tanimoto_similarity_according_to_condition(vec1, vec2, threshold=0.5):
   difference = np.abs(np.array(vec1)-np.array(vec2))
   intersection = 0
   union = 0
   for v1, v2, d in zip(vec1, vec2, difference):
      if d < threshold:
         intersection += 1
      else:
         union += 1
   union += intersection
   return intersection / union

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

#Tanimoto coefficient calculations based on topological indices for connected_graphs_with_7_vertices 
#and exporting the obtained data as pkl and xlsx files.
topological_indices_connected_graphs_with_7_vertices=[]
for i in range(0,len(graphs_with_7_vertices)-1):
   for j in range(i+1,len(graphs_with_7_vertices)):
      vector1=[]
      vector2=[]
      degree_based_indices1=Degree_based_indices(graphs_with_7_vertices[i])
      centrality_measures1=centrality_measures(graphs_with_7_vertices[i])
      energy_and_spectral_radius1=energy_and_spectral_radius_of_graph(graphs_with_7_vertices[i])
      vector1.append(round(degree_based_indices1[0],5))
      vector1.append(round(degree_based_indices1[1],5))
      vector1.append(round(degree_based_indices1[2],5))
      vector1.append(round(degree_based_indices1[3],5))
      vector1.append(round(degree_based_indices1[4],5))
      vector1.append(round(degree_based_indices1[5],5))
      vector1.append(round(degree_based_indices1[6],5))
      vector1.append(round(degree_based_indices1[7],5))
      vector1.append(round(degree_based_indices1[8],5))
      vector1.append(round(Mostar_Szeged_index(graphs_with_7_vertices[i])[0],5))
      vector1.append(round(Mostar_Szeged_index(graphs_with_7_vertices[i])[1],5))
      vector1.append(round(variance_of_graph(graphs_with_7_vertices[i]),5))
      vector1.append(round(total_irregularity(graphs_with_7_vertices[i]),5))
      vector1.append(round(Wiener_and_Balaban_indices(graphs_with_7_vertices[i])[0],5))
      vector1.append(round(Wiener_and_Balaban_indices(graphs_with_7_vertices[i])[1],5))
      vector1.append(round(energy_and_spectral_radius1[0],5))
      vector1.append(round(energy_and_spectral_radius1[1],5))
      vector1.append(round(centrality_measures1[0],5))
      vector1.append(round(centrality_measures1[1],5))
      vector1.append(round(centrality_measures1[2],5))

      
      degree_based_indices2=Degree_based_indices(graphs_with_7_vertices[j])
      centrality_measures2=centrality_measures(graphs_with_7_vertices[j])
      energy_and_spectral_radius2=energy_and_spectral_radius_of_graph(graphs_with_7_vertices[j])
      vector2.append(round(degree_based_indices2[0],5))
      vector2.append(round(degree_based_indices2[1],5))
      vector2.append(round(degree_based_indices2[2],5))
      vector2.append(round(degree_based_indices2[3],5))
      vector2.append(round(degree_based_indices2[4],5))
      vector2.append(round(degree_based_indices2[5],5))
      vector2.append(round(degree_based_indices2[6],5))
      vector2.append(round(degree_based_indices2[7],5))
      vector2.append(round(degree_based_indices2[8],5))
      vector2.append(round(Mostar_Szeged_index(graphs_with_7_vertices[j])[0],5))
      vector2.append(round(Mostar_Szeged_index(graphs_with_7_vertices[j])[1],5))
      vector2.append(round(variance_of_graph(graphs_with_7_vertices[j]),5))
      vector2.append(round(total_irregularity(graphs_with_7_vertices[j]),5))
      vector2.append(round(Wiener_and_Balaban_indices(graphs_with_7_vertices[j])[0],5))
      vector2.append(round(Wiener_and_Balaban_indices(graphs_with_7_vertices[j])[1],5))
      vector2.append(round(energy_and_spectral_radius2[0],5))
      vector2.append(round(energy_and_spectral_radius2[1],5))
      vector2.append(round(centrality_measures2[0],5))
      vector2.append(round(centrality_measures2[1],5))
      vector2.append(round(centrality_measures2[2],5))
      topological_indices_connected_graphs_with_7_vertices.append(tanimoto_similarity_according_to_condition(vector1,vector2))    

df2=pd.DataFrame(topological_indices_connected_graphs_with_7_vertices)
df2.to_excel('topological_indices_connected_graphs_with_7_vertices.xlsx',index=False)

with open("topological_indices_connected_graphs_with_7_vertices.pkl", "wb") as topological_indicesfile:
    pickle.dump(topological_indices_connected_graphs_with_7_vertices, topological_indicesfile)
#