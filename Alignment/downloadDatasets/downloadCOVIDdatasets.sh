#!/bin/bash

module load sratoolkit-precomp
module load ANACONDA/5.0.1
source /apps/ANACONDA/5.0.1/bin/activate rnavel

declare -a samples=("SRR11181954" "SRR11181955" "SRR11181956" "SRR11181957" "SRR11181958" "SRR11181959" "SRR11537946" "SRR11537947" "SRR11537948" "SRR11537949" "SRR11537950" "SRR115379451")



for i in "${samples[@]}"; do
   if [[ ! -d $i ]]; then
     mkdir $i
   fi

   files=$(ls | grep "$i" | wc -l)
   while [[ ! $files -eq 2 ]]; do 
     echo "Downloading sample ${i}"
     prefetch $i --max-size 8000000000 -v -O /gpfs/scratch/bsc08/bsc08991/RNAseqProtocol/code
     parallel-fastq-dump --sra-id $i --threads 48 --split-files --gzip
   done

   echo "Sample ${i} downloaded"
   mv $i_1.fastq.gz ../inputFiles/$1/
   mv $i_2.fastq.gz ../inputFiles/$1/
done
