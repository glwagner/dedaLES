#!/bin/bash
### Job Name
#PBS -N test
### Project code
#PBS -A UMIT0023
#PBS -l walltime=01:00:00
#PBS -q regular
### Merge output and error files
#PBS -j oe
### Select 2 nodes with 36 CPUs each for a total of 72 MPI processes
#PBS -l select=2:ncpus=36:mpiprocs=36
### Send email on abort, begin and end
#PBS -m abe
### Specify mail recipient
#PBS -M wagner.greg@gmail.com

export TMPDIR=/glade/scratch/$USER/temp
mkdir -p $TMPDIR

dedales="$HOME/dedaLES"
examples="$dedales/examples"
scriptname="$examples/rayleigh_benard_example.py"

### Activate Miniconda
. /glade/u/home/$USER/software/miniconda3/etc/profile.d/conda.sh

conda activate dedalus

module load intel/17.0.1 
module load openmpi/3.0.1

### Run the executable
mpiexec python3 $scriptname >> test.out
