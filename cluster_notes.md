# Notes on running `dedaLES` on high-performance computing clusters

The first thing to note is that `dedaLES` has high memory requirements. 
The scale of computational jobs is limited by memory requirements.

## Notes on Engaging

The main partition on MIT's engaging cluster is `sched_mit_hill`.
It has the properties:

```
          Nodes : 228
  CPUs per node : 16
Memory per node : 64000 MB
 Memory per CPU : 4000 MB
```

To specify the maximum memory per node, include the line

```
#SBATCH --mem-per-cpu=4000
```

in your `slurm` script.
This directive tells `slurm` to permit 4 GB of memory allocation per CPU, which is the 
maximum memory per CPU allowed on `sched_mit_hill`.

An example `slurm` script is

```bash
#!/bin/bash

# Job
#SBATCH --partition=sched_mit_hill
#SBATCH --nodes=24
#SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=4000
#SBATCH --time=4:00:00
#SBATCH --job-name="test"

# Streams
#SBATCH --output=job_%j.out
#SBATCH --error=job_%j.err

# Setup (for example, conda activate dedalus)

# Content
mpiexec python3 program.py >> test.out
```

where `program.py` is a `dedaLES` script.
