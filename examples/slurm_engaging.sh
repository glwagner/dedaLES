#!/bin/bash

# Partition             Nodes   S-C-T   Timelimit
# ---------             -----   -----   ---------
# sched_mit_hill        (32)    2-8-1   12:00:00
# sched_mit_raffaele    (32)    2-10-1  12:00:00
# sched_any_quicktest   2       2-8-1   00:15:00
# newnodes              (32)    2-10-1  12:00:00

# Job
#SBATCH --partition=sched_mit_hill
#SBATCH --nodes=24
#SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=4000
#SBATCH --time=4:00:00
#SBATCH --job-name="FC0"

# Streams
#SBATCH --output=job_%j.out
#SBATCH --error=job_%j.err

# Activate conda and dedalus environment
. /home/glwagner/software/miniconda3/etc/profile.d/conda.sh
conda activate dedalus

# Content
mpiexec python3 free_convection_example.py >> FC0.out

#analysis="freeconvection_nx128_ny16_nz128_Q0p00000080544088_Ninv10000_DNS"
#mpiexec python3 merge.py $analysis
