#!/bin/bash
#SBATCH -J parallel_optimize_BTSP_CA1
#SBATCH -o parallel_optimize_BTSP_CA1_102317.%j.o
#SBATCH -e parallel_optimize_BTSP_CA1_102317.%j.e
#SBATCH --partition=compute
#SBATCH -N 17 -n 402
#SBATCH -t 12:00:00
##SBATCH --mail-user=graceyng@stanford.edu
#SBATCH --mail-user=aaronmil@stanford.edu
#SBATCH --mail-type=END
#SBATCH --mail-type=BEGIN
#

set -x

cd $HOME/BTSP
export DATE=$(date +%Y%m%d_%H%M%S)
export SCRATCH=/oasis/scratch/comet/$USER/temp_project
cluster_id="BTSP_CA1_$DATE"

ibrun -n 1 ipcontroller --ip='*' --nodb --cluster-id=$cluster_id &
sleep 1
sleep 60
ibrun -n 400 ipengine --mpi=mpi4py --cluster-id=$cluster_id &
sleep 1
sleep 180
ibrun -n 1 python parallel_optimize.py --config-file-path='data/parallel_optimize_BTSP_CA1_config.yaml' --cluster-id=$cluster_id --pop-size=200 --max-iter=50 --path-length=3 --disp --export --output-dir=$SCRATCH
