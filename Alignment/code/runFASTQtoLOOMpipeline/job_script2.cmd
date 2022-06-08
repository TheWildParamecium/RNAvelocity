#!/bin/bash
#SBATCH --job-name="RNA_vel2"
#SBATCH --workdir=.
#SBATCH --output=test_%j.out
#SBATCH --error=test_%j.err
#SBATCH --cpus-per-task=48
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --time=14:00:00

#./downloadDataset.sh
./fromFASTQtoBAM2.sh
