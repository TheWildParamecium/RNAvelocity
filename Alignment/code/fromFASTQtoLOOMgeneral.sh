#!/bin/bash

####################### Loading modules and configuring enviroment

#########Loading required programs from the HPC modules
module load cellranger/6.0.1
module load samtools

########If you want to reproduce the same workflow, you might need to configure
#an anaconda enviroment like that

module load ANACONDA/5.0.1

isRNAvelEnviroment=$(conda env list | grep "rnavel" | wc -l)
if [[ ${isRNAvelEnviroment} -eq 1 ]]; then
  echo "Anaconda enviroment for RNA velocity is already configured"
else
  conda env create -f requirements.yml -n rnavel
fi

source /apps/ANACONDA/5.0.1/bin/activate rnavel

############################# Applying Cellranger count to align fastq reads to reference genome and produce bam files

### Performing some prior steps to check if the input files are in the correct format

file1=$(ls ../inputFiles/$1/$2 | grep -E ".+\_S[0-9]_L[0-9]+_R1_[0-9]+\.fastq\.gz" | wc -l) 
file2=$(ls ../inputFiles/$1/$2 | grep -E ".+\_S[0-9]_L[0-9]+_R2_[0-9]+\.fastq\.gz" | wc -l) 
filess=$(($file1 + $file2))

if [[ ! `expr $filess % 2` == 0 ]]; then
 echo "Your files are not named correctly. You need to have your sample names (R1 and R2 files) like this:"
 echo "[samplename]_S1_L00[lane number]_R1_001.fastq.gz"
 echo "[samplename]_S1_L00[lane number]_R2_001.fastq.gz"
else
 echo "your input filenames are correctly named"
fi

####################### Running CellRanger script to produce BAM files from FASTQ files

mkdir -p $2/outs
bamfiles=$(ls ./$2/outs | grep "bam" | wc -l) 
if [[ bamfiles -ge 3 ]]; then
  echo "It seems you already applied CellRanger to this dataset before, so skipping this step and procceding with velocytoCLI"
else
  rm -rf ./$2
  cellranger count --id=$2 \
        --transcriptome=../refGenomes/refdata-gex-GRCh38-2020-A \
        --fastqs=../inputFiles/$1/$2 \
        --expect-cells=$3 \
        --localcores=48 \
        --localmem=88 \
        --localvmem=32
fi

###################### Applying velocytoCLI to produce loom files (a file with unspliced and spliced count matrices) from BAM files obtained before

if [[ ! -f ./$2/velocyto/$2 ]]; then
	velocyto run10x -m ../refGenomes/hg38_repeat_msk.gtf \
		./$2 \
		../refGenomes/refdata-gex-GRCh38-2020-A/genes/genes.gtf
else
	echo "It seems like you already applied the workflow to this dataset before, so skipping annotation with velocytoCLI"
fi

########################## Finally, we move the files to the output folder
##(actually, the only interesting file for RNA velocity is the [name].loom file inside the ../outputFiles/[project]/sample/velocyto folder)

if [[ ! -d ../outputFiles/$1/$2 ]]; then
  mkdir -p ../outputFiles/$1/$2
fi

mv ./$2 ../outputFiles/$1

if [[ -f ../outputFiles/${1}/${2}/velocyto/${2}.loom ]]; then
 echo "If this line gets printed, all ran smoothly. You can find your loom file inside ../outputFiles/${1}/${2}/velocyto/${2}.loom"
else
 echo "It seems that something went wrong. Look at the log files, located at ./logs/${1}/${2}"
fi
