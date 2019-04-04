#!/bin/bash
##Run the MIXCR program to find immuneSEq from RNAseq in TCGA pancreatic data
##https://mixcr.readthedocs.io/en/latest/rnaseq.html#ref-rna-seq

# Add relevant executables to PATH
export PATH=$PATH:~/scripts/mixcr-2.1.11/ # MIXCR tool
export PATH=$PATH:~/scripts/gdc-client-dir/ # gdc transfer client

# Set the file input
# input="~/UUID_test.txt" # location of text file containing UUIDs of samples, each on a separate line
token="gdc-user-token.2018-10-13T00_35_47.766Z.txt" # controlled access token
#data="~/data/" # location where all MIXCR data (alignments, etc.) will be located
#TMP_FILE_PREFIX="~/tmp/" # location to temporarily store output of MIXCR commands (raw output, not reports)

num=5 # number of samples that will be processed simultaneously at any given time

# Download and process the TCGA samples in batches of n samples each
count=0
dataset="" # for storing name of each data set
# temp="" # for storing .fastq filename prefix ([...]_1.fastq or [...]_2.fastq)

for input in ~/inputs/*
do
  dataset=$(basename "$input")
  dataset=${dataset::-13} # remove file path at beginning and "_gdc_manifest.2018-07-06.txt0X" at end of name
  echo "CURRENT: $dataset"

  mkdir ~/data/${dataset} # create a directory for current dataset
  data=~/data/${dataset} # location where all MIXCR data (alignments, etc.) will be located

  mkdir ~/tmp/${dataset}
  tmp_storage=~/tmp/${dataset}

  # process the dataset
  while read -r UUID
  do
    # UUID=${UUID::-1} # get rid of "\r" since inputs use carriage return at end of line
    echo $UUID

    # 1. Download sample
    gdc-client download $UUID -t ${token}

    # Wait every ${num} samples for processing to complete such that next sample is downloaded
    # while current n samples are being processed
    if [ $count -eq $num ]
    then
      wait
      echo
      count=0
    fi

    ((count++)) # increment count by 1
    mv $UUID ${data} # move into data directory

    for filename in ${data}/$UUID/*.tar.gz # Untar data
    do
      tar -xzf $filename -C ${data}/$UUID/
    done

    echo "downloaded and expanded"


    # 2. Run MIXCR in background
    (echo "processing ${UUID}" && \
    mixcr align -p rna-seq -s hsa -f -r ${data}/${UUID}/${UUID}_alignments_report.txt -t 36 -OallowPartialAlignments=true `find ${data}/${UUID}/ -name *.fastq` ${data}/${UUID}/${UUID}_alignments.vdjca > ${tmp_storage}/MIXCR_${UUID}.txt && \
    echo "initial alignments complete for ${UUID}" && \
    mixcr assemblePartial -f -r ${data}/${UUID}/${UUID}alignments_report_rescued_1.txt ${data}/${UUID}/${UUID}_alignments.vdjca ${data}/${UUID}/${UUID}_alignments_rescued_1.vdjca >> ${tmp_storage}/MIXCR_${UUID}.txt && \
    mixcr assemblePartial -f -r ${data}/${UUID}/${UUID}alignments_report_rescued_2.txt ${data}/${UUID}/${UUID}_alignments_rescued_1.vdjca ${data}/${UUID}/${UUID}_alignments_rescued_2.vdjca >> ${tmp_storage}/MIXCR_${UUID}.txt && \
    echo "partial alignments complete for ${UUID}" && \
    mixcr extendAlignments -f -r ${data}/${UUID}/${UUID}alignments_report_extended.txt  ${data}/${UUID}/${UUID}_alignments_rescued_2.vdjca ${data}/${UUID}/${UUID}_alignments_extended.vdjca >> ${tmp_storage}/MIXCR_${UUID}.txt && \
    echo "alignments extended for ${UUID}" && \
    mixcr assemble -r ${data}/${UUID}/${UUID}clonotypes_report.txt -i ${data}/${UUID}/${UUID}_index_file ${data}/${UUID}/${UUID}_alignments_extended.vdjca ${data}/${UUID}/${UUID}_output.clns >> ${tmp_storage}/MIXCR_${UUID}.txt && \
    echo "clonotypes assembled for ${UUID}" && \
    mixcr exportAlignments -readId -sequence -vHit -dHit -jHit -vGene -dGene -jGene -vAlignment -dAlignment -jAlignment -nFeature VGene -nFeature CDR3 -aaFeature CDR3 -lengthOf VGene ${data}/${UUID}/${UUID}_alignments_extended.vdjca ${data}/${UUID}/${UUID}_alignments.txt >> ${tmp_storage}/MIXCR_${UUID}.txt && \
    echo "alignments exported for ${UUID}" && \
    mixcr exportClones -f -cloneId -sequence -count -vHit -jHit -vAlignment -jAlignment -nFeature CDR3 -aaFeature CDR3 ${data}/${UUID}/${UUID}_output.clns ${data}/${UUID}/${UUID}_clones.txt >> ${tmp_storage}/MIXCR_${UUID}.txt && \
    echo "clones exported for ${UUID}" && \
    rm `find ${data}/${UUID}/ -name *.fastq` &&\
    echo ".fastq files deleted for ${UUID}" && \
    rm `find ${data}/${UUID}/ -name *.tar.gz` && \
    echo ".tar.gz file deleted for ${UUID}" && \
    aws s3 cp ${data}/${UUID} s3://katyu/mixcr_output/${dataset}/ --recursive && \
    echo "${UUID} COMPLETED AND TRANSFERED" && \
    rm -r ${data}/${UUID} && \
    echo "DELETED ${UUID}") &

  done < "$input"

  wait
  echo
  echo "ALL SAMPLES FOR $(basename $dataset) PROCESSED"
  echo
  echo
done

echo
echo "SCRIPT COMPLETE"
