#!/bin/bash

fetch_seq_files(){
   SraAccList_FILE=${DIR}/SraAccList.txt.pilot
   $SRA_PREFETCH_SCRIPT --option-file $SraAccList_FILE

   echo "Log|  Split .sra files to FW/RV .fastq files"
   for i in $(cat $SraAccList_FILE)
   do
       $SRA_FASTQ_DUMP_SCRIPT -I --split-files  --gzip ./$i/$i.sra
   done
}

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
. ${DIR}/../../../settings.txt
echo "Log| Fetching $DIR"
pushd ${DIR}

fetch_seq_files

popd