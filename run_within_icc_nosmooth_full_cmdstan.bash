#!/bin/bash

sbatch -c 16 --mem=10G -t 1-00:00:00 --array="${1}" ~/data/containers/sbatch_R_command.bash \
	verse-cmdstan.simg \
	within_person_icc.R \
	--chains 4 --threads 4 --adaptdelta .9999 --nosmooth
