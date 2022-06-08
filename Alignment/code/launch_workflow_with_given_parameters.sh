#!/bin/bash

while getopts "n:J:c:t:q" arg; do
  case $arg in
   	c) ncores=${OPTARG} ;;
	n) ntasks=${OPTARG} ;;
	t) time=${OPTARG} ;;
	J) jobname=${OPTARG} ;;
	q) debugging='true' ;;

  esac
done

ARG1=${@:$OPTIND:1}
ARG2=${@:$OPTIND+1:1}
ARG3=${@:$OPTIND+2:1}

if [[ ! -d ./logs/${ARG1}/${ARG2} ]]; then mkdir -p ./logs/${ARG1}/${ARG2}; fi

#Launching sbatch script overwritting some parameters fo the slurm configuration file

if [[ $debugging ]]; then
	sbatch	-n $ntasks -c $ncores -t $time -J $jobname --qos="debug" \
		-e ./logs/${ARG1}/${ARG2}/test_%x_%j.err -o ./logs/${ARG1}/${ARG2}/test_%x_%j.out \
		 jobscript_general.cmd ${ARG1} ${ARG2} ${ARG3}
else
	 sbatch	-n $ntasks -c $ncores -t $time -J $jobname \
		-e ./logs/${ARG1}/${ARG2}/test_%x_%j.err -o ./logs/${ARG1}/${ARG2}/test_%x_%j.out \
		jobscript_general.cmd ${ARG1} ${ARG2} ${ARG3}
fi
