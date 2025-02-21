# measuring_similarity
Data for paper on a novel method for measuring similarity or distance for arbitrary graphs

In the ‘data_acquisition’ file, Tanimoto coefficients of all isomer pairs are calculated separately for undecane, dodecane and tridecane isomers.
Among these codes, the code for undecanes uses the files ‘isomer11.pkl’ and ‘isomers_of_undecanes.xlsx’ as input.
Among these codes, the code for dodecanes uses the files ‘isomer12.pkl’ and ‘isomers_of_dodecanes.xlsx’ as input.
Among these codes, the code for tridecanes uses the files ‘isomer13.pkl’ and ‘isomers_of_tridecanes.xlsx’ as input.
Each of these codes provides the values of the Tanimoto coefficients calculated based on each fingerprint and adjaceny matrix as output in the form of xlsx and pkl files, 6 files each. The xlsx files have been used for plotting the graphics mentioned below. If desired, pkl files can also be used.

The codes in the ‘plots’ file use the output obtained as a result of the codes in the ‘data_acquisition’ file for the undecane, dodecane and tridecane isomers to plot the value distribution plots of the Tanimoto coefficients, the cumulative distribution plots of the Tanimoto coefficients and the relationships between the Tanimoto coefficients based on molecular fingerprints and adjacency matrices.
