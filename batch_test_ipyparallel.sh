#!/bin/bash
#SBATCH -N 3 -n 72
#SBATCH -t 00:30:00
#SBATCH -J 071017_test_ipyparallel
#SBATCH -o 071017_test_ipyparallel.%j.o -e 071017_test_ipyparallel.%j.e
#SBATCH --mail-user=aaronmil@stanford.edu
#SBATCH --mail-type=END
#SBATCH --mail-type=BEGIN

cd $HOME/CA1Sim
cluster_id="$1"

set -x
mpirun -n 1 ipcontroller --ip='*' --cluster-id=$cluster_id &
sleep 60
mpirun -n 70 ipengine --cluster-id=$cluster_id &
sleep 60
mpirun -n 1 python test_ipyparallel.py --cluster-id=$cluster_id
