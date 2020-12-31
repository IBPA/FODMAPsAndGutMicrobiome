#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
pushd $DIR

echo "Log| Fetching ./data/shared/"
wget http://www2.decipher.codes/Classification/TrainingSets/SILVA_SSU_r132_March2018.RData
wget https://gg-sg-web.s3-us-west-2.amazonaws.com/downloads/greengenes_database/gg_13_5/gg_13_5.fasta.gz
wget https://gg-sg-web.s3-us-west-2.amazonaws.com/downloads/greengenes_database/gg_13_5/gg_13_5_otus.tar.gz
tar xzf gg_13_5_otus.tar.gz
rm -rf ./gg_13_5_otus/rep_set_aligned/ # big directory which is not needed now.

## Activate Qiime2 envirnment
. ${DIR}/../../settings.txt
. $CONDA_ENV
conda activate $QIIME_ENV

qiime tools import \
    --type 'FeatureData[Sequence]' \
    --input-path $DIR/gg_13_5_otus/rep_set/97_otus.fasta \
    --output-path $DIR/gg_97_otus.qza

qiime tools import \
    --type 'FeatureData[Sequence]' \
    --input-path $DIR/gg_13_5_otus/rep_set/94_otus.fasta \
    --output-path $DIR/gg_94_otus.qza

qiime tools import \
    --type 'FeatureData[Sequence]' \
    --input-path $DIR/gg_13_5_otus/rep_set/91_otus.fasta \
    --output-path $DIR/gg_91_otus.qza

# After picrust installaion:
download_picrust_files.py

popd

conda deactivate
