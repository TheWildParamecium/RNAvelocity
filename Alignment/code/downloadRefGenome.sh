#!/bin/bash


#For Human genome
###############################################################################

#Uncomment to download the reference genome you want

########### Getting human reference genome for aligment for CellRanger
# wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
# tar -zxvf ./refdata-gex-GRCh38-2020-A.tar.gz
# mv refdata-gex-GRCh38-2020-A ../refGenomes/refdata-gex-GRCh38-2020-A

#For Mouse Genome

########### Getting mouse reference genome for aligment for CellRanger
#wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz
#tar -zxvf ./refdata-gex-mm10-2020-A.tar.gz
#mv refdata-gex-mm10-2020-A ../refGenomes/refdata-gex-mm10-2020-A


##############################################################################


#For shared genome (xenografts experiments)

########### Getting mouse reference genome for aligment for CellRanger
#wget wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-and-mm10-2020-A.tar.gz
#tar -zxvf ./refdata-gex-GRCh38-and-mm10-2020-A.tar.gz
mv refdata-gex-GRCh38-and-mm10-2020-A ../refGenomes/refdata-gex-GRCh38-and-mm10-2020-A
