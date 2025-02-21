# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 15:31:28 2025

@author: sinan.oz
"""

import cirpy as cp
from pysmiles import read_smiles
import networkx as nx
import pickle
import pandas as pd
from rdkit.Chem import AllChem
from typing import Set
from rdkit.Chem import rdMolDescriptors
from rdkit.Avalon import pyAvalonTools
from rdkit.Chem import MACCSkeys

with open("./Desktop/obtaining_data/dodecanes/isomer12.pkl", "rb") as file12:
    isomer12= pickle.load(file12)

#Sorting the graphs in the isomer12 file according to the SMILES order in the isomers_of_dodecanes file
df=pd.read_excel('isomers_of_dodecanes.xlsx')
graphs_of_dodecanes=[None] * 355
for i in range(len(df)):
   smi=cp.resolve(df.Isomer[i], 'smiles')
   G = read_smiles(smi)
   for j in range(len(isomer12)):
       is_isomorphic=nx.is_isomorphic(G, isomer12[j])
       if is_isomorphic:
         graphs_of_dodecanes[i]=G

#Tanimoto coefficient calculations based on the adjacency matrix for dodecanes
#and exporting the obtained data as pkl and xlsx files.
adjacency_matrix_dodecanes= [] 
for i in range(0,len(graphs_of_dodecanes)-1):
    for j in range(i+1,len(graphs_of_dodecanes)):
      Matrix1 = nx.adjacency_matrix(graphs_of_dodecanes[i])
      M1 = Matrix1.toarray()
      Matrix2 = nx.adjacency_matrix(graphs_of_dodecanes[j])
      M2 = Matrix2.toarray()
      res_common_1=0
      res_1_M1=0
      res_1_M2=0
      for i2 in range(0,len(M1)):
         res_common_1 += sum(x == y==1 for x, y in zip(M1[i2], M2[i2]))
         res_1_M1 += sum(z ==1 for z in M1[i2])
         res_1_M2 += sum(y ==1 for y in M2[i2])
      tanimoto_similarity=res_common_1/(res_1_M1+res_1_M2-res_common_1)
      adjacency_matrix_dodecanes.append(tanimoto_similarity)

df2=pd.DataFrame(adjacency_matrix_dodecanes)
df2.to_excel('adjacency_matrix_dodecanes.xlsx',index=False)

with open("adjacency_matrix_dodecanes.pkl", "wb") as adjacency_matrix_dodecanesfile6:
    pickle.dump(adjacency_matrix_dodecanes, adjacency_matrix_dodecanesfile6)
#

def tanimoto_coefficient(set_1: Set[int], set_2: Set[int]) -> float:
    intersect= len(set_1.intersection(set_2))
    union = len(set_1) + len(set_2) - intersect
    return intersect / union

#Tanimoto coefficient calculations based on atom pair for dodecanes
#and exporting the obtained data as pkl and xlsx files.
rows = 355
cols = 512 
N_atom_pair_dodecanes= [[None for _ in range(cols)] for _ in range(rows)]
array=[]
for i in range(len(df)):
   smi=cp.resolve(df.Isomer[i], 'smiles')
   m1 = AllChem.MolFromSmiles(smi)
   AP = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(m1, nBits=512)
   q1=0
   for j in range(len(AP)):
       N_atom_pair_dodecanes[i][q1]=AP[j]
       q1=q1+1

atom_pair_dodecanes= []     
for i in range(0,354):
    for j in range(i+1,355):
       set_a = set(index for index, value in enumerate(N_atom_pair_dodecanes[i]) if value == 1)
       set_b = set(index for index, value in enumerate(N_atom_pair_dodecanes[j]) if value == 1)
       tanimoto_similarity = tanimoto_coefficient(set_a, set_b)
       atom_pair_dodecanes.append(tanimoto_similarity)

df2=pd.DataFrame(atom_pair_dodecanes)
df2.to_excel('atom_pair_dodecanes.xlsx',index=False)

with open("atom_pair_dodecanes.pkl", "wb") as atom_pair_dodecanesfile:
    pickle.dump(atom_pair_dodecanes, atom_pair_dodecanesfile)
#

#Tanimoto coefficient calculations based on topological torsion for dodecanes
#and exporting the obtained data as pkl and xlsx files.
rows = 355
cols = 512 
N_topological_torsion= [[None for _ in range(cols)] for _ in range(rows)]
array=[]
for i in range(len(df)):
   smi=cp.resolve(df.Isomer[i], 'smiles')
   m1 = AllChem.MolFromSmiles(smi)
   TT = rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(m1, nBits=512)
   q1=0
   for j in range(len(TT)):
       N_topological_torsion[i][q1]=TT[j]
       q1=q1+1

topological_torsion_dodecanes= []
for i in range(0,354):
    for j in range(i+1,355):
       set_x = set(index for index, value in enumerate(N_topological_torsion[i]) if value == 1)
       set_y = set(index for index, value in enumerate(N_topological_torsion[j]) if value == 1)
       tanimoto_similarity = tanimoto_coefficient(set_x, set_y)
       topological_torsion_dodecanes.append(tanimoto_similarity)

df2=pd.DataFrame(topological_torsion_dodecanes)
df2.to_excel('topological_torsion_dodecanes.xlsx',index=False)

with open("topological_torsion_dodecanes.pkl", "wb") as topological_torsion_dodecanesfile:
    pickle.dump(topological_torsion_dodecanes, topological_torsion_dodecanesfile)
#

#Tanimoto coefficient calculations based on Avalon for dodecanes
#and exporting the obtained data as pkl and xlsx files.
rows = 355
cols = 512 
N_Avalon_dodecanes= [[None for _ in range(cols)] for _ in range(rows)]
array=[]
for i in range(len(df)):
   smi=cp.resolve(df.Isomer[i], 'smiles')
   m1 = AllChem.MolFromSmiles(smi)
   Av = pyAvalonTools.GetAvalonFP(m1, nBits=512)
   q1=0
   for j in range(len(Av)):
       N_Avalon_dodecanes[i][q1]=Av[j]
       q1=q1+1

Avalon_dodecanes= []     
for i in range(0,354):
    for j in range(i+1,355):
       set_a = set(index for index, value in enumerate(N_Avalon_dodecanes[i]) if value == 1)
       set_b = set(index for index, value in enumerate(N_Avalon_dodecanes[j]) if value == 1)
       tanimoto_similarity = tanimoto_coefficient(set_a, set_b)
       Avalon_dodecanes.append(tanimoto_similarity)

df2=pd.DataFrame(Avalon_dodecanes)
df2.to_excel('Avalon_dodecanes.xlsx',index=False)

with open("Avalon_dodecanes.pkl", "wb") as Avalon_dodecanesfile:
    pickle.dump(Avalon_dodecanes, Avalon_dodecanesfile)
#
    
#Tanimoto coefficient calculations based on MACCS keys for dodecanes
#and exporting the obtained data as pkl and xlsx files.  
rows = 355
cols = 512 
N_MACCS_keys= [[None for _ in range(cols)] for _ in range(rows)]
array=[]
for i in range(len(df)):
   smi=cp.resolve(df.Isomer[i], 'smiles')
   m1 = AllChem.MolFromSmiles(smi)
   Mkeys = MACCSkeys.GenMACCSKeys(m1)
   q1=0
   for j in range(len(Mkeys)):
       N_MACCS_keys[i][q1]=Mkeys[j]
       q1=q1+1

MACCS_keys_dodecanes= []     
for i in range(0,354):
    for j in range(i+1,355):
       set_a = set(index for index, value in enumerate(N_MACCS_keys[i]) if value == 1)
       set_b = set(index for index, value in enumerate(N_MACCS_keys[j]) if value == 1)
       tanimoto_similarity = tanimoto_coefficient(set_a, set_b)
       MACCS_keys_dodecanes.append(tanimoto_similarity)

df2=pd.DataFrame(MACCS_keys_dodecanes)
df2.to_excel('MACCS_keys_dodecanes.xlsx',index=False)

with open("MACCS_keys_dodecanes.pkl", "wb") as MACCS_keys_dodecanesfile:
    pickle.dump(MACCS_keys_dodecanes, MACCS_keys_dodecanesfile)
#

#Tanimoto coefficient calculations based on Morgan circular for dodecanes
#and exporting the obtained data as pkl and xlsx files.
rows = 355
cols = 512 
N_Morgan_circular= [[None for _ in range(cols)] for _ in range(rows)]
array=[]
for i in range(len(df)):
   smi=cp.resolve(df.Isomer[i], 'smiles')
   m1 = AllChem.MolFromSmiles(smi)
   onbits = {}
   mf = AllChem.GetMorganFingerprintAsBitVect(m1, 2, nBits=512, bitInfo=onbits)
   q1=0
   for j in range(len(mf)):
       N_Morgan_circular[i][q1]=mf[j]
       q1=q1+1

Morgan_circular_dodecanes= []
for i in range(0,354):
    for j in range(i+1,355):
       set_a = set(index for index, value in enumerate(N_Morgan_circular[i]) if value == 1)
       set_b = set(index for index, value in enumerate(N_Morgan_circular[j]) if value == 1)
       tanimoto_similarity = tanimoto_coefficient(set_a, set_b)
       Morgan_circular_dodecanes.append(tanimoto_similarity)

df2=pd.DataFrame(Morgan_circular_dodecanes)
df2.to_excel('Morgan_circular_dodecanes.xlsx',index=False)

with open("Morgan_circular_dodecanes.pkl", "wb") as Morgan_circular_dodecanesfile:
    pickle.dump(Morgan_circular_dodecanes, Morgan_circular_dodecanesfile)
#

      
