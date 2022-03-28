#!/bin/bash
#SBATCH -J balpha
#SBATCH --time=3-00:00:00
#SBATCH -n 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH -p ncf
#SBATCH --account=mclaughlin_lab
# Outputs ----------------------------------
#SBATCH -o log/%x-%A_%a.out
#parcel numbers: 108,196,197,235,236,293,305,307,311,312,344,395,396,397,410,420
#Already run: 108,196,197,311,312,397,410,420
#New: 235,236,293,305,307,344,395,396
source /users/jflournoy/code/R_3.5.1_modules.bash

srun -c 4 Rscript --vanilla within_person_icc.R "$@"
