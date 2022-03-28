#!/bin/bash

sbatch -c 16 --mem=32G -t 2-00:00:00 -w holy7c22310 --array="${1}" ~/data/containers/sbatch_R_command.bash \
	verse-cmdstan.simg \
	within_person_icc.R \
	--subsample --chains 4 --threads 4 --adaptdelta .99999
