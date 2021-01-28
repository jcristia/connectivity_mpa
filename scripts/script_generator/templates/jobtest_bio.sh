#!/bin/bash

#SBATCH --time=0-00:30        # time (DD-HH:MM)
#SBATCH --mem-per-cpu=4000M   # memory; default unit is megabytes

#SBATCH --mail-user=jcristia10@gmail.com
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

#SBATCH --job-name=bio_test_1
#SBATCH --output=./outputlogs/%x.out

module load singularity
singularity exec --home /home/jcristia/scratch/mpaconn biology_mpaconn.sif python scripts/simtest/bio_1.py