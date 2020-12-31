#!/bin/bash

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
. ${SCRIPT_DIR}/../settings.txt

function print_usage(){
    echo "script usage: $(basename $0) [-s shared_dir] [-p perc_identity] data_dir" >&2
}

function print_success(){
    echo "PASSED: $1"
}

function print_failure(){
    echo "FAILED: $1"
}

## Parse options Start ###
while getopts ':s:p:' OPTION; do
  case "$OPTION" in
    s)
      SHARED_DATA_DIR=$OPTARG
      ;;
    p)
      PERC_IDENTITY=$OPTARG
      ;;
    ?)
      print_usage
      exit 1
      ;;
  esac
done
shift $((OPTIND -1))

if [ -z "$SHARED_DATA_DIR" ]
then
    SHARED_DATA_DIR=${SCRIPT_DIR}/shared
fi

if [ ! -d "$SHARED_DATA_DIR" ]
then
    echo "directory $SHARED_DATA_DIR does not exist!" >&2
    print_usage
    exit 1
fi

if [ -z "$PERC_IDENTITY" ]
then
    PERC_IDENTITY=97
fi

DATA_DIR=$1
if [ ! -d "$DATA_DIR" ]
then
    echo "directory $DATA_DIR does not exist!"
    print_usage
    exit 1
fi
### Parse options End ###

## Convert DADA2 related files (produced by a DADA2.R) for PICRUSt
${SCRIPT_DIR}/../lib/convert_dada2_out.R \
    -i $DATA_DIR/V1_seqtab.nochim.rds \
    -b $DATA_DIR/seqtab.biom.tsv \
    -f $DATA_DIR/seqtab.fasta\
    --taxa_in $DATA_DIR/V2_taxa.rds\
    --taxa_out $DATA_DIR/taxa_metadata.txt\
    && { print_success "convert_dada2_out.R"  ; }\
    || { print_failure "convert_dada2_out.R" ; exit 1; }

## Activate Qiime2 envirnment
. $CONDA_ENV
conda activate $QIIME_ENV\
    && { print_success "conda activate $QIIME_ENV"  ; }\
    || { print_failure "conda activate $QIIME_ENV" ; exit 1; }

## Convert
biom convert \
    -i $DATA_DIR/seqtab.biom.tsv \
    -o $DATA_DIR/seqtab.biom --to-hdf5\
    && { print_success "biom convert"  ; }\
    || { print_failure "biom convert" ; exit 1; }
biom add-metadata \
    -i $DATA_DIR/seqtab.biom \
    -o $DATA_DIR/seqtab_tax.biom \
    --observation-metadata-fp  $DATA_DIR/taxa_metadata.txt \
    --observation-header OTUID,taxonomy \
    --sc-separated taxonomy\
    && { print_success "biom add-metadata"  ; }\
    || { print_failure "biom add-metadata" ; exit 1; }

## Summary information:
biom summarize-table \
    -i $DATA_DIR/seqtab.biom \
    -o $DATA_DIR/seqtab_summary.txt\
    && { print_success "biom summarize-table"  ; }\
    || { print_failure "biom summarize-table" ; exit 1; }

## Convert files to the Qiime format (.qza)
qiime tools import \
    --type 'FeatureData[Sequence]'\
    --input-path $DATA_DIR/seqtab.fasta\
    --output-path $DATA_DIR/seqtabSeq.qza\
    && { print_success "qiime tools import (Sequence)"  ; }\
    || { print_failure "qiime tools import (Sequence)" ; exit 1; }
qiime tools import \
    --type 'FeatureTable[Frequency]'\
    --input-path $DATA_DIR/seqtab.biom\
    --output-path $DATA_DIR/seqtabFreq.qza\
    && { print_success "qiime tools import (Frequency)"  ; }\
    || { print_failure "qiime tools import (Frequency)" ; exit 1; }

## Perform OTU picking (closed-reference) using Greengenes ref db.
qiime vsearch cluster-features-closed-reference \
  --i-table $DATA_DIR/seqtabFreq.qza \
  --i-sequences $DATA_DIR/seqtabSeq.qza \
  --i-reference-sequences $SHARED_DATA_DIR/gg_${PERC_IDENTITY}_otus.qza \
  --p-perc-identity 0.${PERC_IDENTITY} \
  --p-strand both \
  --o-clustered-table $DATA_DIR/table-cr-${PERC_IDENTITY}.qza \
  --o-clustered-sequences $DATA_DIR/rep-seqs-cr-${PERC_IDENTITY}.qza \
  --o-unmatched-sequences $DATA_DIR/unmatched-cr-${PERC_IDENTITY}.qza\
  && { print_success "qiime vsearch cluster-features-closed-reference"  ; }\
  || { print_failure "qiime vsearch cluster-features-closed-reference" ; exit 1; }

## Remove samples with no Greengene matching OTUs
qiime feature-table filter-samples \
  --i-table $DATA_DIR/table-cr-${PERC_IDENTITY}.qza \
  --p-min-frequency 1 \
  --o-filtered-table $DATA_DIR/table-cr-${PERC_IDENTITY}_filtered.qza\
  && { print_success "qiime feature-table filter-samples"  ; }\
  || { print_failure "qiime feature-table filter-samples" ; exit 1; }

## Convert Greengene OTU frequencies to .biom format for PICRUSt
qiime tools export \
    --input-path $DATA_DIR/table-cr-${PERC_IDENTITY}_filtered.qza \
    --output-path $DATA_DIR\
    && { print_success "qiime tools export"  ; }\
    || { print_failure "qiime tools export" ; exit 1; }
mv $DATA_DIR/feature-table.biom $DATA_DIR/table-cr-${PERC_IDENTITY}_filtered.biom

conda deactivate

## ************ PICRUSt ************
## Activate PICRUSt environment
conda activate $PICRUSt_ENV\
    && { print_success "conda activate $PICRUSt_ENV"  ; }\
    || { print_failure "conda activate $PICRUSt_ENV" ; exit 1; }

## Normalize OTU copy numbers (bcz some OTUs are easier to detect)
normalize_by_copy_number.py \
        -i $DATA_DIR/table-cr-${PERC_IDENTITY}_filtered.biom \
        -o $DATA_DIR/PICRUSt_table-cr-${PERC_IDENTITY}_filtered_normalized.biom\
        && { print_success "normalize_by_copy_number.py" ; } \
        || { print_failure "normalize_by_copy_number.py" ; exit 1; }

## Predict metagenome using PICRUSt
predict_metagenomes.py \
        -i $DATA_DIR/PICRUSt_table-cr-${PERC_IDENTITY}_filtered_normalized.biom \
        -o $DATA_DIR/PICRUSt_metagenome_predictions.biom \
        -a $DATA_DIR/PICRUSt_nsti_per_sample.tsv \
        && { print_success "predict_metagenomes.py" ; } \
        || { print_failure "predict_metagenomes.py" ; exit 1; }
biom convert \
    -i $DATA_DIR/PICRUSt_metagenome_predictions.biom \
    -o $DATA_DIR/PICRUSt_metagenome_predictions.tsv --to-tsv \
    && { print_success "biom convert" ; } \
    || { print_failure "biom convert" ; exit 1; }

## Extract higher level KEGG_Pathways
categorize_by_function.py -f \
    -i $DATA_DIR/PICRUSt_metagenome_predictions.biom \
    -c KEGG_Pathways -l 3 \
    -o $DATA_DIR/PICRUSt_predicted_metagenomes_KEGG_L3.tsv \
    && { print_success "categorize_by_function.py" ; } \
    || { print_failure "categorize_by_function.py" ; exit 1; }
categorize_by_function.py \
    -i $DATA_DIR/PICRUSt_metagenome_predictions.biom \
    -c KEGG_Pathways -l 3 \
    -o $DATA_DIR/PICRUSt_predicted_metagenomes_KEGG_L3.biom \
    && { print_success "categorize_by_function.py" ; } \
    || { print_failure "categorize_by_function.py" ; exit 1; }
    