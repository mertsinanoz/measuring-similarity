# A method for measuring similarity or distance for arbitrary graphs
Data for paper on a novel method for measuring similarity or distance for arbitrary graphs

# Requirements
The packages required to run the scripts are specified in the ‘requirements.txt’ file.

# Usage and Explanations
The files ‘data_acquisition_undecanes’, ‘data_acquisition_dodecanes’ and ‘data_acquisition_tridecanes’ contain codes to calculate the Jaccard/Tanimoto indices based on five different fingerprints and the Jaccard/Tanimoto indices based on the newly defined topological index vector of all isomer pairs separately for undecane, dodecane and tridecane isomers, respectively.

In the ‘data_acquisition_undecanes’, the code for undecanes uses the files ‘isomer11.pkl’ and ‘isomers_of_undecanes.xlsx’ as input.

In the ‘data_acquisition_dodecanes’, the code for dodecanes uses the files ‘isomer12.pkl’ and ‘isomers_of_dodecanes.xlsx’ as input.

In the ‘data_acquisition_tridecanes’, the code for tridecanes uses the files ‘isomer13.pkl’ and ‘isomers_of_tridecanes.xlsx’ as input.

In order for the inputs of the codes mentioned above to be used in the code, the file path must be written as required on line 21 in the codes.

Each of these codes provides the values of the Jaccard/Tanimoto indices calculated based on each fingerprint and newly defined topological index vector as output in the form of xlsx and pkl files, 6 files each. These outputs obtained as a result of running the codes are also uploaded to the repository. The xlsx files have been used for plotting the graphics mentioned below. If desired, pkl files can also be used.

The codes in the ‘plots_for_...’ files use the output obtained as a result of the codes in the ‘data_acquisition_...’ files for the undecane, dodecane and tridecane isomers to plot the value distribution plots of the Jaccard/Tanimoto indices, the cumulative distribution plots of the Jaccard/Tanimoto indices and the relationships between the Jaccard/Tanimoto indices based on molecular fingerprints and newly defined topological index vectors.

In addition, the file 'Supplementary_materials' presents the order of occurrence of the undecane, dodecane and tridecane isomers used in the codes and the values of the Jaccard/Tanimoto indices generated as a result of the codes.
