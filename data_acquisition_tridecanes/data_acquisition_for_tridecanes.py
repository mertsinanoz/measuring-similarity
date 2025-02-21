# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 15:34:18 2025

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

with open("./Desktop/obtaining_data/tridecanes/isomer13.pkl", "rb") as file13:
    isomer13= pickle.load(file13)

#Sorting the graphs in the isomer13 file according to the SMILES order in the isomers_of_tridecanes file
df=pd.read_excel('isomers_of_tridecanes.xlsx')
graphs_of_tridecanes=[None] * 802
for i in range(len(df)):
   smi=cp.resolve(df.Isomer[i], 'smiles')
   G = read_smiles(smi)
   for j in range(len(isomer13)):
       is_isomorphic=nx.is_isomorphic(G, isomer13[j])
       if is_isomorphic:
         graphs_of_tridecanes[i]=G

#Tanimoto coefficient calculations based on the adjacency matrix for tridecanes
#and exporting the obtained data as pkl and xlsx files.
adjacency_matrix_tridecanes= [] 
for i in range(0,len(graphs_of_tridecanes)-1):
    for j in range(i+1,len(graphs_of_tridecanes)):
      Matrix1 = nx.adjacency_matrix(graphs_of_tridecanes[i])
      M1 = Matrix1.toarray()
      Matrix2 = nx.adjacency_matrix(graphs_of_tridecanes[j])
      M2 = Matrix2.toarray()
      res_common_1=0
      res_1_M1=0
      res_1_M2=0
      for i2 in range(0,len(M1)):
         res_common_1 += sum(x == y==1 for x, y in zip(M1[i2], M2[i2]))
         res_1_M1 += sum(z ==1 for z in M1[i2])
         res_1_M2 += sum(y ==1 for y in M2[i2])
      tanimoto_similarity=res_common_1/(res_1_M1+res_1_M2-res_common_1)
      adjacency_matrix_tridecanes.append(tanimoto_similarity)

df2=pd.DataFrame(adjacency_matrix_tridecanes)
df2.to_excel('adjacency_matrix_tridecanes.xlsx',index=False)

with open("adjacency_matrix_tridecanes.pkl", "wb") as adjacency_matrix_tridecanesfile6:
    pickle.dump(adjacency_matrix_tridecanes, adjacency_matrix_tridecanesfile6)
#

def tanimoto_coefficient(set_1: Set[int], set_2: Set[int]) -> float:
    intersect= len(set_1.intersection(set_2))
    union = len(set_1) + len(set_2) - intersect
    return intersect / union

#Tanimoto coefficient calculations based on atom pair for tridecanes
#and exporting the obtained data as pkl and xlsx files.
rows = 802
cols = 512 
N_atom_pair_tridecanes= [[None for _ in range(cols)] for _ in range(rows)]
array=[]
for i in range(len(df)):
   smi=cp.resolve(df.Isomer[i], 'smiles')
   m1 = AllChem.MolFromSmiles(smi)
   AP = rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(m1, nBits=512)
   q1=0
   for j in range(len(AP)):
       N_atom_pair_tridecanes[i][q1]=AP[j]
       q1=q1+1

atom_pair_tridecanes= []     
for i in range(0,801):
    for j in range(i+1,802):
       set_a = set(index for index, value in enumerate(N_atom_pair_tridecanes[i]) if value == 1)
       set_b = set(index for index, value in enumerate(N_atom_pair_tridecanes[j]) if value == 1)
       tanimoto_similarity = tanimoto_coefficient(set_a, set_b)
       atom_pair_tridecanes.append(tanimoto_similarity)

df2=pd.DataFrame(atom_pair_tridecanes)
df2.to_excel('atom_pair_tridecanes.xlsx',index=False)

with open("atom_pair_tridecanes.pkl", "wb") as atom_pair_tridecanesfile:
    pickle.dump(atom_pair_tridecanes, atom_pair_tridecanesfile)
#

#Tanimoto coefficient calculations based on topological torsion for tridecanes
#and exporting the obtained data as pkl and xlsx files.
rows = 802
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

topological_torsion_tridecanes= []
for i in range(0,801):
    for j in range(i+1,802):
       set_x = set(index for index, value in enumerate(N_topological_torsion[i]) if value == 1)
       set_y = set(index for index, value in enumerate(N_topological_torsion[j]) if value == 1)
       tanimoto_similarity = tanimoto_coefficient(set_x, set_y)
       topological_torsion_tridecanes.append(tanimoto_similarity)

df2=pd.DataFrame(topological_torsion_tridecanes)
df2.to_excel('topological_torsion_tridecanes.xlsx',index=False)

with open("topological_torsion_tridecanes.pkl", "wb") as topological_torsion_tridecanesfile:
    pickle.dump(topological_torsion_tridecanes, topological_torsion_tridecanesfile)
#

#Tanimoto coefficient calculations based on Avalon for tridecanes
#and exporting the obtained data as pkl and xlsx files.
rows = 802
cols = 512 
N_Avalon_tridecanes= [[None for _ in range(cols)] for _ in range(rows)]
array=[]
for i in range(len(df)):
   smi=cp.resolve(df.Isomer[i], 'smiles')
   m1 = AllChem.MolFromSmiles(smi)
   Av = pyAvalonTools.GetAvalonFP(m1, nBits=512)
   q1=0
   for j in range(len(Av)):
       N_Avalon_tridecanes[i][q1]=Av[j]
       q1=q1+1

Avalon_tridecanes= []     
for i in range(0,801):
    for j in range(i+1,802):
       set_a = set(index for index, value in enumerate(N_Avalon_tridecanes[i]) if value == 1)
       set_b = set(index for index, value in enumerate(N_Avalon_tridecanes[j]) if value == 1)
       tanimoto_similarity = tanimoto_coefficient(set_a, set_b)
       Avalon_tridecanes.append(tanimoto_similarity)

df2=pd.DataFrame(Avalon_tridecanes)
df2.to_excel('Avalon_tridecanes.xlsx',index=False)

with open("Avalon_tridecanes.pkl", "wb") as Avalon_tridecanesfile:
    pickle.dump(Avalon_tridecanes, Avalon_tridecanesfile)
#
   
#Tanimoto coefficient calculations based on MACCS keys for tridecanes
#and exporting the obtained data as pkl and xlsx files.   
rows = 802
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

MACCS_keys_tridecanes= []     
for i in range(0,801):
    for j in range(i+1,802):
       set_a = set(index for index, value in enumerate(N_MACCS_keys[i]) if value == 1)
       set_b = set(index for index, value in enumerate(N_MACCS_keys[j]) if value == 1)
       tanimoto_similarity = tanimoto_coefficient(set_a, set_b)
       MACCS_keys_tridecanes.append(tanimoto_similarity)

df2=pd.DataFrame(MACCS_keys_tridecanes)
df2.to_excel('MACCS_keys_tridecanes.xlsx',index=False)

with open("MACCS_keys_tridecanes.pkl", "wb") as MACCS_keys_tridecanesfile:
    pickle.dump(MACCS_keys_tridecanes, MACCS_keys_tridecanesfile)
#

#Tanimoto coefficient calculations based on Morgan circular for tridecanes
#and exporting the obtained data as pkl and xlsx files.
rows = 802
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

Morgan_circular_tridecanes= []
for i in range(0,801):
    for j in range(i+1,802):
       set_a = set(index for index, value in enumerate(N_Morgan_circular[i]) if value == 1)
       set_b = set(index for index, value in enumerate(N_Morgan_circular[j]) if value == 1)
       tanimoto_similarity = tanimoto_coefficient(set_a, set_b)
       Morgan_circular_tridecanes.append(tanimoto_similarity)

df2=pd.DataFrame(Morgan_circular_tridecanes)
df2.to_excel('Morgan_circular_tridecanes.xlsx',index=False)

with open("Morgan_circular_tridecanes.pkl", "wb") as Morgan_circular_tridecanesfile:
    pickle.dump(Morgan_circular_tridecanes, Morgan_circular_tridecanesfile)
#