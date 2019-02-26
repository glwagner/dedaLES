#!/bin/bash

# Partition             Nodes   S-C-T   Timelimit
# ---------             -----   -----   ---------
# sched_mit_hill        (32)    2-8-1   12:00:00
# sched_mit_raffaele    (32)    2-10-1  12:00:00
# sched_any_quicktest   2       2-8-1   00:15:00
# newnodes              (32)    2-10-1  12:00:00

# Job
#SBATCH --partition=sched_mit_hill
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=3500
#SBATCH --time=12:00:00
#SBATCH --job-name="kerr1-CS"

# Streams
#SBATCH --output=job_%j.out
#SBATCH --error=job_%j.err

# Activate conda and dedalus environment
. /home/glwagner/software/miniconda3/etc/profile.d/conda.sh
conda activate dedalus

run="1"
closure="ConstantSmagorinsky"

# Content
mpiexec python3 rayleigh_benard_kerr.py $run $closure >> kerr_CS_$run.out
