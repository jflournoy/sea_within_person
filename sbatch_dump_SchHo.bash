#!/bin/bash
#SBATCH -J dump
#SBATCH --time=0-01:00:00
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3G
#SBATCH -p ncf
#SBATCH --account=mclaughlin_lab
# Outputs ----------------------------------
#SBATCH -o log/%x-%A_%a.out

module load afni/18.3.05-fasrc01

set -aeuxo pipefail

N=$SLURM_ARRAY_TASK_ID

#/ContrastName1  Calm > Fear 
#/ContrastName2  Calm > Happy 
#/ContrastName3  Calm > Scramble 
#/ContrastName4  Fear > Calm 
#/ContrastName5  Fear > Happy 
#/ContrastName6  Fear > Scramble 
#/ContrastName7  Happy > Calm 
#/ContrastName8  Happy > Fear 
#/ContrastName9  Happy > Scramble 
#/ContrastName10 Calm 
#/ContrastName11 Fear 
#/ContrastName12 Happy 
#/ContrastName13 Scramble 

# 292 /users/jflournoy/NewNeuropoint/datafiles/setfilenames_fearGTcalm.txt
files=($(cat /mnt/stressdevlab/stress_pipeline/Group/FaceReactivity/NewNeuropoint/datafiles/setfilenames_fearGTcalm.txt))
file=${files[$N]}

calm=${file%cope4.nii.gz}cope10.nii.gz
fear=${file%cope4.nii.gz}cope11.nii.gz
happy=${file%cope4.nii.gz}cope12.nii.gz
happyGTcalm=${file%cope4.nii.gz}cope7.nii.gz
copes=(${file} ${calm} ${fear} ${happy} ${happyGTcalm})
names=(fearGTcalm calm fear happy happyGTcalm)

outfilepre=$(echo $file | sed -e 's/.*\([0-9]\{4\}\)\/month\([0-9]\{2\}\).*/\1_\2/')

J=(0 1 2 3 4)
for jj in ${J[@]}; do
  if [ ! -f  "tsv/${outfilepre}_${names[${jj}]}.tsv" ]; then
    3dmaskdump -mask Sch400_HO_combined.nii.gz -nozero -o "tsv/${outfilepre}_${names[${jj}]}.tsv" "${copes[${jj}]}"
    fi
done;
