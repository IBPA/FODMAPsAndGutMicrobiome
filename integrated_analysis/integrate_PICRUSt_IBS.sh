#!/bin/bash

function print_success(){
    echo "PASSED: $1"
}

function print_failure(){
    echo "FAILED: $1"
}

function convert_metagenome_KEGG(){
    DATA_DIR=$1
    qiime tools import \
        --type 'FeatureTable[Frequency]'\
        --input-path $DATA_DIR/PICRUSt_predicted_metagenomes_KEGG_L3.biom\
        --output-path $DATA_DIR/PICRUSt_predicted_metagenomes_KEGG_L3.qza\
        && { print_success "qiime tools import for $DATA_DIR"  ; }\
        || { print_failure "qiime tools import for $DATA_DIR" ; exit 1; }
}

function convert_metagenome_all(){
    DATA_DIR=$1
    qiime tools import \
        --type 'FeatureTable[Frequency]'\
        --input-path $DATA_DIR/PICRUSt_metagenome_predictions.biom\
        --output-path $DATA_DIR/PICRUSt_metagenome_predictions.qza\
        && { print_success "qiime tools import for $DATA_DIR"  ; }\
        || { print_failure "qiime tools import for $DATA_DIR" ; exit 1; }
}

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
. ${DIR}/../settings.txt
. $CONDA_ENV
conda activate $QIIME_ENV

OUT_DIR=${DIR}/integrated_IBS
mkdir $OUT_DIR

#### KEGG Metagenome ####
S1_micr_data_dir=${DIR}/../preprocessing/s1/microbiome_data/
S2_micr_data_dir=${DIR}/../preprocessing/s2/microbiome_data/
S3_micr_data_dir=${DIR}/../preprocessing/s3/microbiome_data/

convert_metagenome_KEGG ${S1_micr_data_dir}
convert_metagenome_KEGG ${S2_micr_data_dir}
convert_metagenome_KEGG ${S3_micr_data_dir}

qiime feature-table merge\
    --i-tables ${S1_micr_data_dir}/PICRUSt_predicted_metagenomes_KEGG_L3.qza\
    --i-tables ${S2_micr_data_dir}/PICRUSt_predicted_metagenomes_KEGG_L3.qza\
    --i-tables ${S3_micr_data_dir}/PICRUSt_predicted_metagenomes_KEGG_L3.qza\
    --o-merged-table ${OUT_DIR}/PICRUSt_predicted_metagenomes_KEGG_L3.qza\
    && { print_success "qiime feature-table merge"  ; }\
    || { print_failure "qiime feature-table merge" ; exit 1; }

qiime tools export \
    --input-path ${OUT_DIR}/PICRUSt_predicted_metagenomes_KEGG_L3.qza \
    --output-path ${OUT_DIR}\
    && { print_success "qiime tools export"  ; }\
    || { print_failure "qiime tools export" ; exit 1; }
mv ${OUT_DIR}/feature-table.biom  ${OUT_DIR}/PICRUSt_predicted_metagenomes_KEGG_L3_NoTaxa.biom

biom add-metadata \
    -i ${OUT_DIR}/PICRUSt_predicted_metagenomes_KEGG_L3_NoTaxa.biom \
    -o ${OUT_DIR}/PICRUSt_predicted_metagenomes_KEGG_L3.biom \
    --observation-metadata-fp  ${DIR}/shared/PICRUSt_KEGG_Pathways.tsv \
    --observation-header OTUID,KEGG_Pathways \
    --sc-separated KEGG_Pathways\
    && { print_success "biom add-metadata"  ; }\
    || { print_failure "biom add-metadata" ; exit 1; }


#### All Metagenome ####

convert_metagenome_all ${S1_micr_data_dir}
convert_metagenome_all ${S2_micr_data_dir}
convert_metagenome_all ${S3_micr_data_dir}

qiime feature-table merge\
    --i-tables ${S1_micr_data_dir}/PICRUSt_metagenome_predictions.qza\
    --i-tables ${S2_micr_data_dir}/PICRUSt_metagenome_predictions.qza\
    --i-tables ${S3_micr_data_dir}/PICRUSt_metagenome_predictions.qza\
    --o-merged-table ${OUT_DIR}/PICRUSt_metagenome_predictions.qza\
    && { print_success "qiime feature-table merge"  ; }\
    || { print_failure "qiime feature-table merge" ; exit 1; }

qiime tools export \
    --input-path ${OUT_DIR}/PICRUSt_metagenome_predictions.qza \
    --output-path ${OUT_DIR}\
    && { print_success "qiime tools export"  ; }\
    || { print_failure "qiime tools export" ; exit 1; }
mv ${OUT_DIR}/feature-table.biom  ${OUT_DIR}/PICRUSt_metagenome_predictions_NoTaxa.biom

biom add-metadata \
    -i ${OUT_DIR}/PICRUSt_metagenome_predictions_NoTaxa.biom \
    -o ${OUT_DIR}/PICRUSt_metagenome_predictions.biom \
    --observation-metadata-fp  ${DIR}/shared/PICRUSt_metagenome_KO.tsv \
    --observation-header OTUID,KEGG_Description,KEGG_Pathways \
    --sc-separated KEGG_Pathways \
    && { print_success "biom add-metadata"  ; }\
    || { print_failure "biom add-metadata" ; exit 1; }

biom convert \
    -i ${OUT_DIR}/PICRUSt_metagenome_predictions.biom \
    -o ${OUT_DIR}/PICRUSt_metagenome_predictions.tsv --to-tsv \
     && { print_success "biom convert"  ; }\
    || { print_failure "biom convert" ; exit 1; }

conda deactivate
