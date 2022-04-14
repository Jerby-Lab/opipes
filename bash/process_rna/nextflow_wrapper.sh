#!/bin/bash

# See `man sbatch` or https://slurm.schedmd.com/sbatch.html for descriptions
# of sbatch options.
#SBATCH --job-name=nextflowwrapper
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=interactive
#SBATCH --account=default
#SBATCH --time=50:00:00

# nextflow/20.01.0
module load nextflow/20.01.0

FASTQ=$1/*.fastq # First argument is the path to the fastq files
EMAIL=$2@stanford.edu # Second argument provided is your SUNetID
echo "Fastq files: $FASTQ"
echo "Email $EMAIL"

# [nextflow] run the rna-seq single end rna! 
nextflow run nf-core/rnaseq -profile singularity \
	--reads $FASTQ \
		--genome hg38 -c rnaseq.config --outdir rna-out --email $email \
			-resume --singleEnd --skipBiotypeQC -w rna-workdir


# How to run? sbatch nextflow_wrapper.sh path_to_fastq SUNetID
# This code is part of the Jerby lab computational resources and documentation: https://github.com/Jerby-Lab/
