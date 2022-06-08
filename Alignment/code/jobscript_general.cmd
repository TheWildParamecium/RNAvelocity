#!/bin/bash
#SBATCH --workdir=.
#SBATCH --job-name="RNAvelocity"
#SBATCH --output=test_%x_%j.out
#SBATCH --error=test_%x_%j.err
#SBATCH --cpus-per-task=48
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --time=48:00:00

./fromFASTQtoLOOMgeneral.sh $1 $2 $3
