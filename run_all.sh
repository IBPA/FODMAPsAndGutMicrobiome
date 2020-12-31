#!/bin/bash

# 1) Perform microbiome data preprocessing for each study
./preprocessing/shared/fetch.sh
./preprocessing/s1/run.sh
./preprocessing/s2/run.sh
./preprocessing/s3/run.sh
./preprocessing/s4/run.sh
./preprocessing/s5/run.sh
./preprocessing/s6/run.sh

# 2) Integrate profiles for IBS studies (16S)
mkdir ./integrated_analysis/integrated_IBS
./integrated_analysis/integrate_PICRUSt_IBS.sh
./integrated_analysis/integrate_DADA2_IBS.R

# 3) Integrate profiles for IBS studies (GA-map)
mkdir ./integrated_analysis/integrated_GA_map
./integrated_analysis/integrate_GA_map.R

# 4) Integrate profiles for all 16S together
mkdir ./integrated_analysis/integrated_all
./integrated_analysis/integrate_PICRUSt_all.sh
./integrated_analysis/integrate_DADA2_all.R

# 5) Generate main figures
mkdir results
./integrated_analysis/generate_Fig2.R
./integrated_analysis/generate_Fig3.R
./integrated_analysis/generate_Fig4.R
./integrated_analysis/generate_Fig5.R

# 6) Generate supplementary figures
./integrated_analysis/generate_FigS1.R
./integrated_analysis/generate_FigS2.R
./integrated_analysis/generate_FigS3.R
./integrated_analysis/generate_FigS4.R
./integrated_analysis/generate_FigS5.R
./integrated_analysis/generate_FigS6.R
./integrated_analysis/generate_FigS78.R
./integrated_analysis/generate_FigS9.R

# 7) Generate supplementary data
./integrated_analysis/generate_som_data.R