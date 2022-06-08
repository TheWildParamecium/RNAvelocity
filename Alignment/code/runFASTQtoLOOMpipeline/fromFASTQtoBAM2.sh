#!/bin/bash


####################### Loading modules and configuring enviroment

#Loading required programs from the HPC modules
module load cellranger/6.0.1
module load samtools

#If you want to reproduce the same workflow, you might need to configure
#an anaconda enviroment like that

#module load anaconda
#conda env create -f requeriments.yml -n rnavel
#conda activate rnavel

########################## Applying workflow for Patient 1 --- 0 days of treament

#cellranger count --id=sampleRNAvel1-0days \
#                   --transcriptome=../refGenomes/refdata-gex-GRCh38-2020-A \
#                   --fastqs=../inputFiles/CLL1-0days \
#                   --expect-cells=3000 \
#                   --localcores=48 \
#                   --localmem=70 \
#		   --localvmem=32

#velocyto run10x -m ../refGenomes/hg38_repeat_msk.gtf ./sampleRNAvel1-0days ../refGenomes/refdata-gex-GRCh38-2020-A/genes/genes.gtf

#mv ./sampleRNAvel1-0days ../outputFiles/sampleRNAvel1-0days

#echo "IF THIS LINE GETS PRINTED IT MEANS THAT ALL THE WORKFLOW WORKED WITHOUT ANY PROBLEMS"

########################## Applying workflow for Patient 1 --- 120 days of treament

#cellranger count --id=sampleRNAvel1-120days \
#                   --transcriptome=../refGenomes/refdata-gex-GRCh38-2020-A \
#                   --fastqs=../inputFiles/CLL1-120days \
#                   --expect-cells=3000 \
#                   --localcores=48 \
#                   --localmem=96 \
#                  --localvmem=32

#velocyto run10x -m ../refGenomes/hg38_repeat_msk.gtf ./sampleRNAvel1-120days ../refGenomes/refdata-gex-GRCh38-2020-A/genes/genes.gtf

#mv ./sampleRNAvel1-120days ../outputFiles/sampleRNAvel1-120days

#echo "IF THIS LINE GETS PRINTED IT MEANS THAT ALL THE WORKFLOW WORKED WITHOUT ANY PROBLEMS"

########################## Applying workflow for human prostate cancer cells in a mouse

#cellranger count --id=mouseprostateRNAvel \
#                   --transcriptome=../refGenomes/refdata-gex-GRCh38-2020-A \
#                   --fastqs=../inputFiles/mouse-prostate \
#                   --expect-cells=3000 \
#                   --localcores=48 \
#                   --localmem=96 \
#                   --localvmem=32 \
#		    --r1-length=26

#velocyto run10x -m ../refGenomes/hg38_repeat_msk.gtf ../outputFiles/mouseprostateRNAvel ../refGenomes/refdata-gex-GRCh38-2020-A/genes/genes.gtf



#mv ./mouseprostateRNAvel ../outputFiles/mouseprostateRNAvel 

#echo "IF THIS LINE GETS PRINTED IT MEANS THAT ALL THE WORKFLOW WORKED WITHOUT ANY PROBLEMS"


############################# Sperm dataset

cellranger count --id=azoospermia34 \
                   --transcriptome=../refGenomes/refdata-gex-GRCh38-2020-A \
                   --fastqs=../inputFiles/azoospermia34 \
                   --expect-cells=6000 \
                   --localcores=48 \
                   --localmem=96 \
                   --localvmem=32 

module load anaconda
conda env create -f requeriments.yml -n rnavel
source /apps/ANACONDA/5.0.1/bin/activate rnavel

velocyto run10x -m ../refGenomes/hg38_repeat_msk.gtf ./azoospermia34 ../refGenomes/refdata-gex-GRCh38-2020-A/genes/genes.gtf
