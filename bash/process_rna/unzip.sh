#!/bin/bash
#SBATCH -J unzip_rna
#SBATCH -t 1:0:0
#SBATCH -o /labs/ccurtis2/cyyeh/tonic/workflows/process_rna/slurm/slurm-"%j".out

# Account 
#SBATCH --account=ccurtis2

# Send emails
#SBATCH --mail-user=cyyeh@stanford.edu
#SBATCH --mail-type=FAIL

# GLOBALS 
WORKDIR='/labs/ccurtis2/cyyeh/tonic/data/RNA-Seq/'

# WORKFLOW
cd $WORKDIR
files=($(ls */*fastq.gz))
file=${files[${SLURM_ARRAY_TASK_ID}-1]}
fastq="${file%.*}"
gunzip -c $file > /labs/ccurtis2/cyyeh/tonic/data/RNA-Seq/$fastq
echo $file
echo $fastq


