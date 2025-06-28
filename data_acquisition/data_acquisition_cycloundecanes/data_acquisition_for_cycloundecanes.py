# -*- coding: utf-8 -*-
"""
Created on Mon Jun  9 16:44:47 2025

@author: sinan.oz
"""

import cirpy as cp
import numpy as np
from pysmiles import read_smiles
import networkx as nx
import pickle
import pandas as pd
import math
from rdkit.Chem import AllChem
from typing import Set
from rdkit.Chem import rdMolDescriptors
from rdkit.Avalon import pyAvalonTools
from rdkit.Chem import MACCSkeys

from rdkit import Chem


def graph_to_smiles(G):
    mol = Chem.RWMol()
    node_to_atom = {}
    for node in G.nodes():
        atom = Chem.Atom('C')
        atom_index = mol.AddAtom(atom)
        node_to_atom[node] = atom_index
    for u, v in G.edges():
        mol.AddBond(node_to_atom[u], node_to_atom[v], Chem.BondType.SINGLE)
    smiles = Chem.MolToSmiles(mol)
    return smiles

unigraphsS11=[]   
with open('./data_acquisition_cycloundecanes/eleven_vertices.g6','r') as f:
     for line in f:
       unigraphsS11.append(line.strip())
       
unicyclicgraphs_with_eleven_vertices=[]     
for i in range(len(unigraphsS11)):
   str1=unigraphsS11[i]
   g1=nx.from_graph6_bytes(bytes(str1, 'ascii', errors='strict'))
   unicyclicgraphs_with_eleven_vertices.append(g1)
    
k = 0
graphs_of_cycloundecanes= []

for tr in range(len(unicyclicgraphs_with_eleven_vertices)):
    M1 = nx.adjacency_matrix(unicyclicgraphs_with_eleven_vertices[tr])
    M = M1.toarray()
    degrees = M.sum(axis=1)
    if all(degrees <= 4):
        k += 1
        graphs_of_cycloundecanes.append(unicyclicgraphs_with_eleven_vertices[tr])
        
smiles_list=[]
for i in range(0,len(graphs_of_cycloundecanes)):
   smi = graph_to_smiles(graphs_of_cycloundecanes[i])
   smiles_list.append(smi)
       
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

#Tanimoto coefficient calculations based on topological indices for cycloundecanes 
#and exporting the obtained data as pkl and xlsx files.
topological_indices_cycloundecanes=[]
for i in range(0,len(graphs_of_cycloundecanes)-1):
   for j in range(i+1,len(graphs_of_cycloundecanes)):
      vector1=[]
      vector2=[]
      degree_based_indices1=Degree_based_indices(graphs_of_cycloundecanes[i])
      centrality_measures1=centrality_measures(graphs_of_cycloundecanes[i])
      energy_and_spectral_radius1=energy_and_spectral_radius_of_graph(graphs_of_cycloundecanes[i])
      vector1.append(round(degree_based_indices1[0],5))
      vector1.append(round(degree_based_indices1[1],5))
      vector1.append(round(degree_based_indices1[2],5))
      vector1.append(round(degree_based_indices1[3],5))
      vector1.append(round(degree_based_indices1[4],5))
      vector1.append(round(degree_based_indices1[5],5))
      vector1.append(round(degree_based_indices1[6],5))
      vector1.append(round(degree_based_indices1[7],5))
      vector1.append(round(degree_based_indices1[8],5))
      vector1.append(round(Mostar_Szeged_index(graphs_of_cycloundecanes[i])[0],5))
      vector1.append(round(Mostar_Szeged_index(graphs_of_cycloundecanes[i])[1],5))
      vector1.append(round(variance_of_graph(graphs_of_cycloundecanes[i]),5))
      vector1.append(round(total_irregularity(graphs_of_cycloundecanes[i]),5))
      vector1.append(round(Wiener_and_Balaban_indices(graphs_of_cycloundecanes[i])[0],5))
      vector1.append(round(Wiener_and_Balaban_indices(graphs_of_cycloundecanes[i])[1],5))
      vector1.append(round(energy_and_spectral_radius1[0],5))
      vector1.append(round(energy_and_spectral_radius1[1],5))
      vector1.append(round(centrality_measures1[0],5))
      vector1.append(round(centrality_measures1[1],5))
      vector1.append(round(centrality_measures1[2],5))

      
      degree_based_indices2=Degree_based_indices(graphs_of_cycloundecanes[j])
      centrality_measures2=centrality_measures(graphs_of_cycloundecanes[j])
      energy_and_spectral_radius2=energy_and_spectral_radius_of_graph(graphs_of_cycloundecanes[j])
      vector2.append(round(degree_based_indices2[0],5))
      vector2.append(round(degree_based_indices2[1],5))
      vector2.append(round(degree_based_indices2[2],5))
      vector2.append(round(degree_based_indices2[3],5))
      vector2.append(round(degree_based_indices2[4],5))
      vector2.append(round(degree_based_indices2[5],5))
      vector2.append(round(degree_based_indices2[6],5))
      vector2.append(round(degree_based_indices2[7],5))
      vector2.append(round(degree_based_indices2[8],5))
      vector2.append(round(Mostar_Szeged_index(graphs_of_cycloundecanes[j])[0],5))
      vector2.append(round(Mostar_Szeged_index(graphs_of_cycloundecanes[j])[1],5))
      vector2.append(round(variance_of_graph(graphs_of_cycloundecanes[j]),5))
      vector2.append(round(total_irregularity(graphs_of_cycloundecanes[j]),5))
      vector2.append(round(Wiener_and_Balaban_indices(graphs_of_cycloundecanes[j])[0],5))
      vector2.append(round(Wiener_and_Balaban_indices(graphs_of_cycloundecanes[j])[1],5))
      vector2.append(round(energy_and_spectral_radius2[0],5))
      vector2.append(round(energy_and_spectral_radius2[1],5))
      vector2.append(round(centrality_measures2[0],5))
      vector2.append(round(centrality_measures2[1],5))
      vector2.append(round(centrality_measures2[2],5))
      topological_indices_cycloundecanes.append(tanimoto_similarity_according_to_condition(vector1,vector2))    

df2=pd.DataFrame(topological_indices_cycloundecanes)
df2.to_excel('topological_indices_cycloundecanes.xlsx',index=False)

with open("topological_indices_cycloundecanes.pkl", "wb") as topological_indicesfile:
    pickle.dump(topological_indices_cycloundecanes, topological_indicesfile)
#

#Tanimoto coefficient calculations based on atom pair for cycloundecanes
#and exporting the obtained data as pkl and xlsx files.
rows = 1231
cols = 512 
N_atom_pair_cycloundecanes= [[None for _ in range(cols)] for _ in range(rows)]
array=[]
for i in range(len(smiles_list)):
   smi=smiles_list[i]
   m1 = AllChem.MolFromSmiles(smi)
   AP = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(m1, nBits=512)
   q1=0
   for j in range(len(AP)):
       N_atom_pair_cycloundecanes[i][q1]=AP[j]
       q1=q1+1

atom_pair_cycloundecanes= []     
for i in range(0,1230):
    for j in range(i+1,1231):
       set_a = set(index for index, value in enumerate(N_atom_pair_cycloundecanes[i]) if value == 1)
       set_b = set(index for index, value in enumerate(N_atom_pair_cycloundecanes[j]) if value == 1)
       tanimoto_similarity = tanimoto_coefficient(set_a, set_b)
       atom_pair_cycloundecanes.append(tanimoto_similarity)

df2=pd.DataFrame(atom_pair_cycloundecanes)
df2.to_excel('atom_pair_cycloundecanes.xlsx',index=False)

with open("atom_pair_cycloundecanes.pkl", "wb") as atom_pair_cycloundecanesfile:
    pickle.dump(atom_pair_cycloundecanes, atom_pair_cycloundecanesfile)
#

#Tanimoto coefficient calculations based on topological torsion for cycloundecanes
#and exporting the obtained data as pkl and xlsx files.
rows = 1231
cols = 512 
N_topological_torsion= [[None for _ in range(cols)] for _ in range(rows)]
array=[]
for i in range(len(smiles_list)):
   smi=smiles_list[i]
   m1 = AllChem.MolFromSmiles(smi)
   TT = rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(m1, nBits=512)
   q1=0
   for j in range(len(TT)):
       N_topological_torsion[i][q1]=TT[j]
       q1=q1+1

topological_torsion_cycloundecanes= []
for i in range(0,1230):
    for j in range(i+1,1231):
       set_x = set(index for index, value in enumerate(N_topological_torsion[i]) if value == 1)
       set_y = set(index for index, value in enumerate(N_topological_torsion[j]) if value == 1)
       tanimoto_similarity = tanimoto_coefficient(set_x, set_y)
       topological_torsion_cycloundecanes.append(tanimoto_similarity)

df2=pd.DataFrame(topological_torsion_cycloundecanes)
df2.to_excel('topological_torsion_cycloundecanes.xlsx',index=False)

with open("topological_torsion_cycloundecanes.pkl", "wb") as topological_torsion_cycloundecanesfile:
    pickle.dump(topological_torsion_cycloundecanes, topological_torsion_cycloundecanesfile)
#

#Tanimoto coefficient calculations based on Avalon for cycloundecanes
#and exporting the obtained data as pkl and xlsx files.
rows = 1231
cols = 512 
N_Avalon_cycloundecanes= [[None for _ in range(cols)] for _ in range(rows)]
array=[]
for i in range(len(smiles_list)):
   smi=smiles_list[i]
   m1 = AllChem.MolFromSmiles(smi)
   Av = pyAvalonTools.GetAvalonFP(m1, nBits=512)
   q1=0
   for j in range(len(Av)):
       N_Avalon_cycloundecanes[i][q1]=Av[j]
       q1=q1+1

Avalon_cycloundecanes= []     
for i in range(0,1230):
    for j in range(i+1,1231):
       set_a = set(index for index, value in enumerate(N_Avalon_cycloundecanes[i]) if value == 1)
       set_b = set(index for index, value in enumerate(N_Avalon_cycloundecanes[j]) if value == 1)
       tanimoto_similarity = tanimoto_coefficient(set_a, set_b)
       Avalon_cycloundecanes.append(tanimoto_similarity)

df2=pd.DataFrame(Avalon_cycloundecanes)
df2.to_excel('Avalon_cycloundecanes.xlsx',index=False)

with open("Avalon_cycloundecanes.pkl", "wb") as Avalon_cycloundecanesfile:
    pickle.dump(Avalon_cycloundecanes, Avalon_cycloundecanesfile)
#
     
#Tanimoto coefficient calculations based on MACCS keys for cycloundecanes
#and exporting the obtained data as pkl and xlsx files.
rows = 1231
cols = 512 
N_MACCS_keys= [[None for _ in range(cols)] for _ in range(rows)]
array=[]
for i in range(len(smiles_list)):
   smi=smiles_list[i]
   m1 = AllChem.MolFromSmiles(smi)
   Mkeys = MACCSkeys.GenMACCSKeys(m1)
   q1=0
   for j in range(len(Mkeys)):
       N_MACCS_keys[i][q1]=Mkeys[j]
       q1=q1+1

MACCS_keys_cycloundecanes= []     
for i in range(0,1230):
    for j in range(i+1,1231):
       set_a = set(index for index, value in enumerate(N_MACCS_keys[i]) if value == 1)
       set_b = set(index for index, value in enumerate(N_MACCS_keys[j]) if value == 1)
       tanimoto_similarity = tanimoto_coefficient(set_a, set_b)
       MACCS_keys_cycloundecanes.append(tanimoto_similarity)

df2=pd.DataFrame(MACCS_keys_cycloundecanes)
df2.to_excel('MACCS_keys_cycloundecanes.xlsx',index=False)

with open("MACCS_keys_cycloundecanes.pkl", "wb") as MACCS_keys_cycloundecanesfile:
    pickle.dump(MACCS_keys_cycloundecanes, MACCS_keys_cycloundecanesfile)
#

#Tanimoto coefficient calculations based on Morgan circular for cycloundecanes
#and exporting the obtained data as pkl and xlsx files.
rows = 1231
cols = 512 
N_Morgan_circular= [[None for _ in range(cols)] for _ in range(rows)]
array=[]
for i in range(len(smiles_list)):
   smi=smiles_list[i]
   m1 = AllChem.MolFromSmiles(smi)
   onbits = {}
   mf = AllChem.GetMorganFingerprintAsBitVect(m1, 2, nBits=512, bitInfo=onbits)
   q1=0
   for j in range(len(mf)):
       N_Morgan_circular[i][q1]=mf[j]
       q1=q1+1

Morgan_circular_cycloundecanes= []
for i in range(0,1230):
    for j in range(i+1,1231):
       set_a = set(index for index, value in enumerate(N_Morgan_circular[i]) if value == 1)
       set_b = set(index for index, value in enumerate(N_Morgan_circular[j]) if value == 1)
       tanimoto_similarity = tanimoto_coefficient(set_a, set_b)
       Morgan_circular_cycloundecanes.append(tanimoto_similarity)

df2=pd.DataFrame(Morgan_circular_cycloundecanes)
df2.to_excel('Morgan_circular_cycloundecanes.xlsx',index=False)

with open("Morgan_circular_cycloundecanes.pkl", "wb") as Morgan_circular_cycloundecanesfile:
    pickle.dump(Morgan_circular_cycloundecanes, Morgan_circular_cycloundecanesfile)
#