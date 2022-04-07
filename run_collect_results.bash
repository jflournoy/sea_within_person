#!/bin/bash

sbatch -c 32 --mem=32G -t 1-00:00:00 ~/data/containers/sbatch_R_command.bash \
	verse-cmdstan.simg \
	collect_results.R
