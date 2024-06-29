# Modeling and analysis code for the E. coli extracellular electron transfer project

ElectricColiModelingClean.m contains scripts for flux balance analysis of the strains and conditions used in this project. These scripts require MATLAB, the COBRA toolbox (https://opencobra.github.io/cobratoolbox/stable/index.html), and a convex solver such as Gurobi (https://www.gurobi.com/) to run.

A preprint of this work is available here: https://www.biorxiv.org/content/10.1101/2024.05.30.596743v1.full

# For the transcriptomics

Unzip all the zipped files in the transcriptomics/models folder before proceeding to run the notebooks

Use requrements_deseq.txt to set up a virtual environment ONLY for the DEGs.ipynb notebook. Use requirements_1.txt to set up a virtual environment to run the other notebooks.
Notebook descriptions:
1) DEGs.ipynb contains code to perform differential gene expression analysis for a pair of samples and get enriched treemaps
2) DiMA_and_Treemaps.ipynb contains code to perform differential imodulon activity analysis for a pair of samples and makes enriched treemaps of imodulon functionality
3) fig5.ipynb contains all the code to generate the figures used in figure 5
3) fig6.ipynb contains all the code to generate the figures used in figure 6

