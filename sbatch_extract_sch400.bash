#!/bin/bash
#SBATCH -J 3droi
#SBATCH --time=0-01:00:00
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3G
#SBATCH -p ncf
#SBATCH --account=mclaughlin_lab
# Outputs ----------------------------------
#SBATCH -o log/%x-%A-%a.out

module load afni/18.3.05-fasrc01

set -aeuxo pipefail

DIRLIST=($( cat $1 ))
N=$SLURM_ARRAY_TASK_ID
SUBSESS=${DIRLIST[${N}]}
outfilepre=$(echo ${SUBSESS} | sed -e 's/.*\([0-9]\{4\}\)\/month\([0-9]\{2\}\).*/\1_\2/')

TASK="faceReactivity"
RUN="FaceReactivity1"
FUNC="${SUBSESS}/${TASK}/${RUN}_new_ssmooth.nii.gz"
SCH400=${SUBSESS}/${TASK}/Schaefer2018_400Parcels_7Networks_order_FSLMNI152_2mm_${RUN}.nii.gz
MASK=${SUBSESS}/${TASK}/${RUN}_brainmask.nii.gz
SCH400NEW=${SUBSESS}/${TASK}/Schaefer2018_400Parcels_7Networks_order_FSLMNI152_2mm_${RUN}_masked.nii.gz

if [ ! -f ${FUNC} ]; then
	echo "File not found: ${FUNC}"
	exit 1
fi

3dcalc -a ${SCH400} -b ${MASK} -expr 'a*b' -prefix ${SCH400NEW}

3dROIstats -mask ${SCH400NEW} -numROI 400 -1Dformat "${FUNC}" > sch400/"${outfilepre}_TS.1D"
