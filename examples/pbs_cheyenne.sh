#!/bin/bash
### Job Name
#PBS -N FC2
### Project code
#PBS -A UMIT0023
#PBS -l walltime=04:00:00
#PBS -q regular
### Merge output and error files
#PBS -j oe
### Select 2 nodes with 36 CPUs each for a total of 72 MPI processes
#PBS -l select=4:ncpus=36:mpiprocs=36:mem=109GB
### Send email on abort, begin and end
#PBS -m abe
### Specify mail recipient
#PBS -M wagner.greg@gmail.com

export TMPDIR=/glade/scratch/$USER/temp
mkdir -p $TMPDIR

export DEDALES="$HOME/dedaLES"

examples="$DEDALES/examples"
benchmark="$DEDALES/benchmarks/rayleigh_benard"
scriptname="free_convection_example.py"

### Activate Miniconda
. /glade/u/home/$USER/software/miniconda3/etc/profile.d/conda.sh

conda activate dedalus

module load intel/17.0.1 
module load openmpi/3.0.1

cp $examples/$scriptname $TMPDIR/$scriptname
cp $examples/merge.py $TMPDIR/
cd $TMPDIR

### Run the executable
mpiexec python3 $scriptname >> $examples/FC2.out

### analysis="freeconvection_nh64_nz64_10Q1_Ninv20000_DNS"
### mpiexec python3 merge.py $analysis
