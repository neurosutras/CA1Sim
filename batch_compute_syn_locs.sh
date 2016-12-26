#!/bin/bash
#
#SBATCH -J neurotrees_full
#SBATCH -o ./data/compute_GC_synapse_locs_122516.%j.o
#SBATCH -n 512
#SBATCH -t 24:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --mail-user=aaronmil@stanford.edu
#SBATCH --mail-type=END
#SBATCH --mail-type=BEGIN
#

set -x
mpirun python compute_GC_synapse_locs.py


