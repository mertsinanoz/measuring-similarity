# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 14:11:44 2025

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

with open("./Desktop/GITHUB/data_acquisition/undecanes/isomer11.pkl", "rb") as file11:
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
#

#Tanimoto coefficient calculations based on the adjacency matrix for undecanes
#and exporting the obtained data as pkl and xlsx files.
adjacency_matrix_undecanes= [] 
for i in range(0,len(graphs_of_undecanes)-1):
    for j in range(i+1,len(graphs_of_undecanes)):
      Matrix1 = nx.adjacency_matrix(graphs_of_undecanes[i])
      M1 = Matrix1.toarray()
      Matrix2 = nx.adjacency_matrix(graphs_of_undecanes[j])
      M2 = Matrix2.toarray()
      res_common_1=0
      res_1_M1=0
      res_1_M2=0
      for i2 in range(0,len(M1)):
         res_common_1 += sum(x == y==1 for x, y in zip(M1[i2], M2[i2]))
         res_1_M1 += sum(z ==1 for z in M1[i2])
         res_1_M2 += sum(y ==1 for y in M2[i2])
      tanimoto_similarity=res_common_1/(res_1_M1+res_1_M2-res_common_1)
      adjacency_matrix_undecanes.append(tanimoto_similarity)

df2=pd.DataFrame(adjacency_matrix_undecanes)
df2.to_excel('adjacency_matrix_undecanes.xlsx',index=False)

with open("adjacency_matrix_undecanes.pkl", "wb") as adjacency_matrix_undecanesfile6:
    pickle.dump(adjacency_matrix_undecanes, adjacency_matrix_undecanesfile6)
#

def tanimoto_coefficient(set_1: Set[int], set_2: Set[int]) -> float:
    intersect= len(set_1.intersection(set_2))
    union = len(set_1) + len(set_2) - intersect
    return intersect / union

#Tanimoto coefficient calculations based on atom pair for undecanes
#and exporting the obtained data as pkl and xlsx files.
rows = 159
cols = 512 
N_atom_pair_undecanes= [[None for _ in range(cols)] for _ in range(rows)]
array=[]
for i in range(len(df)):
   smi=cp.resolve(df.Isomer[i], 'smiles')
   m1 = AllChem.MolFromSmiles(smi)
   AP = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(m1, nBits=512)
   q1=0
   for j in range(len(AP)):
       N_atom_pair_undecanes[i][q1]=AP[j]
       q1=q1+1

atom_pair_undecanes= []     
for i in range(0,158):
    for j in range(i+1,159):
       set_a = set(index for index, value in enumerate(N_atom_pair_undecanes[i]) if value == 1)
       set_b = set(index for index, value in enumerate(N_atom_pair_undecanes[j]) if value == 1)
       tanimoto_similarity = tanimoto_coefficient(set_a, set_b)
       atom_pair_undecanes.append(tanimoto_similarity)

df2=pd.DataFrame(atom_pair_undecanes)
df2.to_excel('atom_pair_undecanes.xlsx',index=False)

with open("atom_pair_undecanes.pkl", "wb") as atom_pair_undecanesfile:
    pickle.dump(atom_pair_undecanes, atom_pair_undecanesfile)
#

#Tanimoto coefficient calculations based on topological torsion for undecanes
#and exporting the obtained data as pkl and xlsx files.
rows = 159
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

topological_torsion_undecanes= []
for i in range(0,158):
    for j in range(i+1,159):
       set_x = set(index for index, value in enumerate(N_topological_torsion[i]) if value == 1)
       set_y = set(index for index, value in enumerate(N_topological_torsion[j]) if value == 1)
       tanimoto_similarity = tanimoto_coefficient(set_x, set_y)
       topological_torsion_undecanes.append(tanimoto_similarity)

df2=pd.DataFrame(topological_torsion_undecanes)
df2.to_excel('topological_torsion_undecanes.xlsx',index=False)

with open("topological_torsion_undecanes.pkl", "wb") as topological_torsion_undecanesfile:
    pickle.dump(topological_torsion_undecanes, topological_torsion_undecanesfile)
#

#Tanimoto coefficient calculations based on Avalon for undecanes
#and exporting the obtained data as pkl and xlsx files.
rows = 159
cols = 512 
N_Avalon_undecanes= [[None for _ in range(cols)] for _ in range(rows)]
array=[]
for i in range(len(df)):
   smi=cp.resolve(df.Isomer[i], 'smiles')
   m1 = AllChem.MolFromSmiles(smi)
   Av = pyAvalonTools.GetAvalonFP(m1, nBits=512)
   q1=0
   for j in range(len(Av)):
       N_Avalon_undecanes[i][q1]=Av[j]
       q1=q1+1

Avalon_undecanes= []     
for i in range(0,158):
    for j in range(i+1,159):
       set_a = set(index for index, value in enumerate(N_Avalon_undecanes[i]) if value == 1)
       set_b = set(index for index, value in enumerate(N_Avalon_undecanes[j]) if value == 1)
       tanimoto_similarity = tanimoto_coefficient(set_a, set_b)
       Avalon_undecanes.append(tanimoto_similarity)

df2=pd.DataFrame(Avalon_undecanes)
df2.to_excel('Avalon_undecanes.xlsx',index=False)

with open("Avalon_undecanes.pkl", "wb") as Avalon_undecanesfile:
    pickle.dump(Avalon_undecanes, Avalon_undecanesfile)
#
     
#Tanimoto coefficient calculations based on MACCS keys for undecanes
#and exporting the obtained data as pkl and xlsx files.
rows = 159
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

MACCS_keys_undecanes= []     
for i in range(0,158):
    for j in range(i+1,159):
       set_a = set(index for index, value in enumerate(N_MACCS_keys[i]) if value == 1)
       set_b = set(index for index, value in enumerate(N_MACCS_keys[j]) if value == 1)
       tanimoto_similarity = tanimoto_coefficient(set_a, set_b)
       MACCS_keys_undecanes.append(tanimoto_similarity)

df2=pd.DataFrame(MACCS_keys_undecanes)
df2.to_excel('MACCS_keys_undecanes.xlsx',index=False)

with open("MACCS_keys_undecanes.pkl", "wb") as MACCS_keys_undecanesfile:
    pickle.dump(MACCS_keys_undecanes, MACCS_keys_undecanesfile)
#

#Tanimoto coefficient calculations based on Morgan circular for undecanes
#and exporting the obtained data as pkl and xlsx files.
rows = 159
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

Morgan_circular_undecanes= []
for i in range(0,158):
    for j in range(i+1,159):
       set_a = set(index for index, value in enumerate(N_Morgan_circular[i]) if value == 1)
       set_b = set(index for index, value in enumerate(N_Morgan_circular[j]) if value == 1)
       tanimoto_similarity = tanimoto_coefficient(set_a, set_b)
       Morgan_circular_undecanes.append(tanimoto_similarity)

df2=pd.DataFrame(Morgan_circular_undecanes)
df2.to_excel('Morgan_circular_undecanes.xlsx',index=False)

with open("Morgan_circular_undecanes.pkl", "wb") as Morgan_circular_undecanesfile:
    pickle.dump(Morgan_circular_undecanes, Morgan_circular_undecanesfile)
#

      
