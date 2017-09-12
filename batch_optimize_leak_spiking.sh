#!/bin/bash
#SBATCH -J parallel_optimize_leak_spiking
#SBATCH -o parallel_optimize_leak_spiking_090417.%j.o
#SBATCH -e parallel_optimize_leak_spiking_090417.%j.e
#SBATCH --partition=compute
#SBATCH -N 42 -n 1002
#SBATCH -t 24:00:00
##SBATCH --mail-user=graceyng@stanford.edu
#SBATCH --mail-user=aaronmil@stanford.edu
#SBATCH --mail-type=END
#SBATCH --mail-type=BEGIN
#

set -x

cd $HOME/CA1Sim
export DATE=$(date +%Y%m%d_%H%M%S)
export SCRATCH=/oasis/scratch/comet/$USER/temp_project
cluster_id="leak_spiking_$DATE"

ibrun -n 1 ipcontroller --ip='*' --DictDB.record_limit=4096 --DictDB.size_limit=3221225472 --DictDB.cull_fraction=0.05 --cluster-id=$cluster_id &
sleep 1
sleep 60
ibrun -n 1000 ipengine --mpi=mpi4py --cluster-id=$cluster_id &
sleep 1
sleep 180
ibrun -n 1 python parallel_optimize.py --param-file-path='data/parallel_optimize_leak_spiking_config.yaml' --cluster-id=$cluster_id --pop-size=200 --max-iter=50 --path-length=3 --disp --export --output-dir=$SCRATCH
