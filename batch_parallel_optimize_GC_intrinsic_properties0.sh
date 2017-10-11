#!/bin/bash
#SBATCH -J parallel_optimize_GC_intrinsic_properties0
#SBATCH -o parallel_optimize_GC_intrinsic_properties0_101117.%j.o
#SBATCH -e parallel_optimize_GC_intrinsic_properties0_101117.%j.e
#SBATCH --partition=compute
#SBATCH -N 42 -n 1002
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
cluster_id="GC_intrinsic0_$DATE"

ibrun -n 1 ipcontroller --ip='*' --nodb --cluster-id=$cluster_id &
sleep 1
sleep 60
ibrun -n 1000 ipengine --mpi=mpi4py --cluster-id=$cluster_id &
sleep 1
sleep 180
ibrun -n 1 python parallel_optimize.py --config-file-path='data/parallel_optimize_GC_intrinsic_properties_config0.yaml' --cluster-id=$cluster_id --pop-size=200 --max-iter=50 --path-length=3 --disp --export --output-dir=$SCRATCH
