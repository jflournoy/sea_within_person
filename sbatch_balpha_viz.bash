#!/bin/bash
#SBATCH -J balpha
#SBATCH --time=3-00:00:00
#SBATCH -n 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -p ncf
#SBATCH --account=mclaughlin_lab
# Outputs ----------------------------------
#SBATCH -o log/%x-%A_%a.out
#parcel numbers: 108,196,197,235,236,293,305,307,311,312,344,395,396,397,410,420
#Already run: 108,196,197,311,312,397,410,420
#New: 235,236,293,305,307,344,395,396
source /users/jflournoy/code/R_3.5.1_modules.bash
#module load rstudio/1.1.453-fasrc01
module load Pandoc/2.5 

srun -c 4 Rscript --vanilla -e "rmarkdown::render('viz.Rmd', output_file = \"report_ROI${SLURM_ARRAY_TASK_ID}_PANDOCTEST.html\", params = list(ROI = ${SLURM_ARRAY_TASK_ID}, SUBSAMP = TRUE) )"
