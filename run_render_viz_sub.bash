#!/bin/bash

sbatch -c 1 --mem=7G -t 1:00:00 --array="${1}" -J "viz" ~/data/containers/sbatch_R_command.bash \
	verse-cmdstan.simg \
	render_viz.R \
	--subsample
