# Biosecurity-based-breeding-schemes
Author: Ziyi Kang kangziyi1998@163.com; Sheng Luan luansheng@ysfri.ac.cn

This code is associated with an unpublished paper titled "Exploring re-ranking magnitude and genomic selection strategies to enhance genetic gain in selective breeding program for Pacific white shrimp (Litopenaeus vannamei) in the presence of genotype by environment interaction".

# Prerequisites
Before running this program, make sure you have the following software and R packages installed on your system:

## Software:

GCTA
Plink
BLUPF90

## R Packages:

AlphaSimR
data.table
visPedigree
optiSel
nadiv
sampling
AlphaMME
AlphaLearn

Please ensure that all the scripts are placed in the same working directory.

# Getting Started
Follow these steps to use this program:

## Set Breeding Scheme Parameters:

Open the parameters.R file and configure the parameters for your breeding schemes. This file controls the specifics of your breeding program.
Set Magnitude of Genotype by Environment Interaction (GEI):

Use the globescript.R script to set the magnitude of Genotype by Environment Interaction. This step allows you to customize the GEI component according to your study's requirements.

## Run the Program:

Execute the globescript.R to run the program. This script will use the parameters set in parameters.R and the GEI settings to perform the calculations and simulations related to your breeding schemes.

# Breeding schemes

bs.R <- BS with genomic selection

ds.R <- NBS with genomic selection

tbs.R <- BS with pedigree-based selection

tds.R <- NBS with pedigree-based selection

# Note
Two scripts, ocs.R and ocsped.R, are sourced from the code of AlphaLearn. These scripts have been included in this project to aid in specific functionalities. Please refer to the original AlphaLearn documentation or source code for additional information regarding these scripts.

For any questions, issues, or inquiries related to this program, feel free to contact [kangziyi1998@163.com].
