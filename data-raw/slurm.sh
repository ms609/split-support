#!/bin/bash

# Request resources:
#SBATCH -c 8            # number of CPU cores, one per thread, up to 128
#SBATCH --mem=24M       # memory required, up to 250G on standard nodes
#SBATCH --gres=tmp:40M  # temporary disk space required on the compute node ($TMPDIR), up to 400G
#SBATCH --time=0:20:0   # time limit for job (format:  days-hours:minutes:seconds)

# Run in the 'shared' queue (job may share node with other jobs)
#SBATCH -p shared

# Commands to be run:
module load bioinformatics
module load MrBayes
