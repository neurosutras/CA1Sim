#!/bin/bash
#SBATCH -J parallel_optimize_GC_intrinsic_properties_scaling1
#SBATCH -o parallel_optimize_GC_intrinsic_properties_scaling1_101117.%j.o
#SBATCH -e parallel_optimize_GC_intrinsic_properties_scaling1_101117.%j.e
#SBATCH --partition=compute
#SBATCH -N 1 -n 12
#SBATCH -t 36:00:00
##SBATCH --mail-user=graceyng@stanford.edu
#SBATCH --mail-user=aaronmil@stanford.edu
#SBATCH --mail-type=END
#SBATCH --mail-type=BEGIN
#

set -x

cd $HOME/CA1Sim
export DATE=$(date +%Y%m%d_%H%M%S)
export SCRATCH=/oasis/scratch/comet/$USER/temp_project
cluster_id="GC_intrinsic_scaling1_$DATE"

ibrun -n 1 ipcontroller --ip='*' --nodb --cluster-id=$cluster_id &
sleep 1
sleep 60
ibrun -n 10 ipengine --mpi=mpi4py --cluster-id=$cluster_id &
sleep 1
sleep 180
ibrun -n 1 python parallel_optimize.py --config-file-path='data/parallel_optimize_GC_intrinsic_properties_config_scaling1.yaml' --cluster-id=$cluster_id --pop-size=10 --max-iter=1 --path-length=1 --disp output-dir=$SCRATCH
