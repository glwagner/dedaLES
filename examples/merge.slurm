#!/bin/bash

#SBATCH --partition=sched_mit_hill
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=4000
#SBATCH --time=1:00:00
#SBATCH --job-name=merge

# Streams
#SBATCH --output=job_%j.out
#SBATCH --error=job_%j.err

# Activate conda and dedalus environment
. /home/glwagner/software/miniconda3/etc/profile.d/conda.sh
conda activate dedalus

# Content
mpiexec python3 merge_script.py freeconvection_nx256_ny256_nz256_F0p0000000001_Ninv1000_DNS >> merge_freeconvection_nx256_ny256_nz256_F0p0000000001_Ninv1000_DNS.out
#!/bin/bash

#SBATCH --reservation=wagner
#SBATCH --partition=sched_mit_hill
#SBATCH --qos=plenum
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=4000
#SBATCH --time=1:00:00
#SBATCH --job-name=merge

# Streams
#SBATCH --output=merge_job_%j.out
#SBATCH --error=merge_job_%j.err

# Activate conda and dedalus environment
. /home/glwagner/software/miniconda3/etc/profile.d/conda.sh
conda activate dedalus

# Content
mpiexec python3 merge_script.py freeconvection_nx256_ny256_nz256_F0p0000000001_Ninv1000_DNS/ >> merge_freeconvection_nx256_ny256_nz256_F0p0000000001_Ninv1000_DNS/.out
#!/bin/bash

#SBATCH --reservation=wagner
#SBATCH --partition=sched_mit_hill
#SBATCH --qos=plenum
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=4000
#SBATCH --time=1:00:00
#SBATCH --job-name=merge

# Streams
#SBATCH --output=merge_job_%j.out
#SBATCH --error=merge_job_%j.err

# Activate conda and dedalus environment
. /home/glwagner/software/miniconda3/etc/profile.d/conda.sh
conda activate dedalus

# Content
mpiexec python3 merge_script.py freeconvection_nx256_ny256_nz256_F0p0000000001_Ninv1000_DNS/ >> merge_freeconvection_nx256_ny256_nz256_F0p0000000001_Ninv1000_DNS/.out
