#!/bin/bash

# Partition             Nodes   S-C-T   Timelimit
# ---------             -----   -----   ---------
# sched_mit_hill        (32)    2-8-1   12:00:00
# sched_mit_raffaele    (32)    2-10-1  12:00:00
# sched_any_quicktest   2       2-8-1   00:15:00
# newnodes              (32)    2-10-1  12:00:00

# Job
#SBATCH --partition=sched_mit_hill
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=3500
#SBATCH --time=2:00:00
###SBATCH --job-name="FC2"
#SBATCH --job-name="merge"

# Streams
#SBATCH --output=job_%j.out
#SBATCH --error=job_%j.err

# Activate conda and dedalus environment
. /home/glwagner/software/miniconda3/etc/profile.d/conda.sh
conda activate dedalus

# Content
mpiexec python3 free_convection_example.py >> FC2.out

analysis="freeconvection_nh128_nz128_10Q10_Ninv10000_DNS"
mpiexec python3 merge.py $analysis
