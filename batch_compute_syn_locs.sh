#!/bin/bash
#
#SBATCH -J neurotrees_full
#SBATCH -o ./results/compute_GC_synapse_locs.%j.o
#SBATCH -n 192
#SBATCH -t 8:00:00
#SBATCH --mem-per-cpu=6G
#SBATCH --mail-user=ivan.g.raikov@gmail.com
#SBATCH --mail-type=END
#SBATCH --mail-type=BEGIN
#
set -x
export LD_LIBRARY_PATH=$PI_HOME/hdf5/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$PI_HOME/python_modules/btmorph:$PYTHONPATH
export SYN_START_INDEX=$1

mpirun python compute_GC_synapse_locs.py


