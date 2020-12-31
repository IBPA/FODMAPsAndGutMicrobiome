#!/bin/bash
study_dir=`dirname "$0"`
pushd $study_dir
./microbiome_data/fetch.sh
./DADA2.R
../PICRUSt_one_study.sh -p 91 $study_dir/microbiome_data
popd