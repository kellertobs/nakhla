#!/bin/bash

#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=120:00:00
#SBATCH --mem=8G

cd ~/nakhla/usr
module load matlab/r2023b

matlab -nosplash -nodisplay -singleCompThread -r "run_$RUNID" > $JOBID.log
