#!/bin/bash

sbatch -c 16 --mem=32G -t 2-00:00:00 --array="${1}" -J "iccFULL" ~/data/containers/sbatch_R_command.bash \
	verse-cmdstan.simg \
	within_person_icc.R \
	--chains 4 --threads 4 
