#!/bin/bash

# Run the full meta-analysis

#######################################
# Run the preprocessing pipeline for a given study
#######################################
function run_preprocessing(){
  study_dir=$1
  echo "Preprocessing study under" $study_dir
  pushd $study_dir
  # ToDo: call the corresponding scripts
  #...
  popd
}

#######################################
# Perform integrated analysis
#######################################
function run_integrated_analysis(){
  echo "Performing integrated analysis and generating the figures."
  ./integrated_analysis/generate_Fig2.R
  ./integrated_analysis/generate_Fig3.R
  ./integrated_analysis/generate_Fig4.R
  ./integrated_analysis/generate_Fig5.R
}

run_preprocessing "./preprocessing/s1"
run_preprocessing "./preprocessing/s2"
run_preprocessing "./preprocessing/s3"
run_preprocessing "./preprocessing/s4"
run_preprocessing "./preprocessing/s5"
run_preprocessing "./preprocessing/s6"

run_integrated_analysis
