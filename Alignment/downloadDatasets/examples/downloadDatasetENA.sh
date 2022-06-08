#!/bin/bash

#################### Patient One --- 0 Days of treatment

#Downloading the files from the ENA ftp repository
# wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR676/005/SRR6762975/SRR6762975_1.fastq.gz
# wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR676/005/SRR6762975/SRR6762975_2.fastq.gz

#Unzipping the files
# gunzip ./SRR6762975_1.fastq.gz
# gunzip ./SRR6762975_2.fastq.gz

#Renaming files to the format CellRanger needs to work
# mv ./SRR6762975_1.fastq ./SRR6762975_S1_L001_R1_001.fastq
# mv ./SRR6762975_2.fastq ./SRR6762975_S1_L001_R2_001.fastq

#Moving the files to the inputFiles directory
# mv ./SRR6762975_S1_L001_R1_001.fastq ../inputFiles/CLL1-0days/SRR6762975_S1_L001_R1_001.fastq
# mv ./SRR6762975_S1_L001_R2_001.fastq ../inputFiles/CLL1-0days/SRR6762975_S1_L001_R2_001.fastq

################## Patient One --- 120 Days after treatment

#Downloading the files from the ENA ftp repository
#wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR676/006/SRR6762976/SRR6762976_1.fastq.gz
#wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR676/006/SRR6762976/SRR6762976_2.fastq.gz

#Unzipping the files
#gunzip ./SRR6762976_1.fastq.gz
#gunzip ./SRR6762976_2.fastq.gz

#Renaming files to the format CellRanger needs to work
# mv ./SRR6762976_1.fastq ./SRR6762976_S1_L001_R1_001.fastq
# mv ./SRR6762976_2.fastq ./SRR6762976_S1_L001_R2_001.fastq

#Moving the files to the inputFiles directory
#mv ./SRR6762976_S1_L001_R1_001.fastq ../inputFiles/CLL1-120days/SRR6762976_S1_L001_R1_001.fastq
#mv ./SRR6762976_S1_L001_R2_001.fastq ../inputFiles/CLL1-120days/SRR6762976_S1_L001_R2_001.fastq



################## (MOUSE) Advanced Prostate Cancer 

#Downloading the files from the ENA ftp repository
#wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR127/068/SRR12700968/SRR12700968_1.fastq.gz
#wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR127/068/SRR12700968/SRR12700968_2.fastq.gz

#Unzipping the files
#gunzip SRR12700968_1.fastq.gz
#gunzip SRR12700968_2.fastq.gz

#Renaming files to the format CellRanger needs to work
mv ./SRR12700968_1.fastq ./SRR12700968_S1_L001_R1_001.fastq
mv ./SRR12700968_2.fastq ./SRR12700968_S1_L001_R2_001.fastq

#Moving the files to the inputFiles directory
#mv ./SRR12700968_S1_L001_R1_001.fastq ../inputFiles/mouse-prostate/SRR12700968_S1_L001_R1_001.fastq
#mv ./SRR12700968_S1_L001_R2_001.fastq ../inputFiles/mouse-prostate/SRR12700968_S1_L001_R2_001.fastq


####################### (to confirm RNA velocity results)
