#!/bin/bash
#SBATCH --workdir=.
#SBATCH --job-name="debug"
#SBATCH --output=test_%x_%j.out
#SBATCH --error=test_%x_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --time=00:10:00
#SBATCH --qos="debug"

./fromFASTQtoLOOMgeneral.sh $1 $2 $3
