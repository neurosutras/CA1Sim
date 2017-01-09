#!/bin/bash
#
#SBATCH -J compute_GC_syn_locs_test
#SBATCH -o ./data/compute_GC_synapse_locs_test_010617.%j.o
#SBATCH -n 512
#SBATCH -t 48:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=1
#BATCH --mail-user=aaronmil@stanford.edu
#BATCH --mail-type=END
#BATCH --mail-type=BEGIN
#
set -x
export LD_LIBRARY_PATH=$PI_HOME/hdf5/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$PI_HOME/python_modules/btmorph:$PYTHONPATH
#export SYN_START_INDEX=$1

mpirun python compute_GC_synapse_locs.py
