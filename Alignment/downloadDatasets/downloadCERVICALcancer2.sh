#!/bin/bash

filesSample1=$(ls ./SRR13927092 | wc -l)
if [[ ! $filesSample1 -eq 2 ]]; then
  mkdir SRR13927092
  wget -c --tries=0 --read-timeout=20 ftp.sra.ebi.ac.uk/vol1/fastq/SRR139/092/SRR13927092/SRR13927092_1.fastq.gz
  wget -c --tries=0 --read-timeout=20 ftp.sra.ebi.ac.uk/vol1/fastq/SRR139/092/SRR13927092/SRR13927092_2.fastq.gz
  mv SRR13927092_1.fastq.gz ./SRR13927092
  mv SRR13927092_2.fastq.gz ./SRR13927092
fi

filesSample2=$(ls ./SRR13927093 | wc -l)
if [[ ! $filesSample2 -eq 2 ]]; then
  mkdir SRR13927093
  wget -c --tries=0 --read-timeout=20 ftp.sra.ebi.ac.uk/vol1/fastq/SRR139/093/SRR13927093/SRR13927093_1.fastq.gz
  wget -c --tries=0 --read-timeout=20 ftp.sra.ebi.ac.uk/vol1/fastq/SRR139/093/SRR13927093/SRR13927093_2.fastq.gz
  mv SRR13927093_1.fastq.gz ./SRR13927093
  mv SRR13927093_2.fastq.gz ./SRR13927093
fi
