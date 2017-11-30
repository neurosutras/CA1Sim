#!/bin/bash -l

#SBATCH -J parallel_optimize_BTSP_CA1
#SBATCH -o parallel_optimize_BTSP_CA1_112917_full.%j.o
#SBATCH -e parallel_optimize_BTSP_CA1_112917_full.%j.e
#SBATCH -p regular
#SBATCH -N 240
#SBATCH -L SCRATCH
#SBATCH -C haswell
#SBATCH -t 18:00:00
#SBATCH --mail-user=aaronmil@stanford.edu
#SBATCH --mail-type=BEGIN,END,FAIL

source $HOME/.bash_profile
set -x

cd $HOME/BTSP
export DATE=$(date +%Y%m%d_%H%M%S)
cluster_id="BTSP_CA1_$DATE"

srun -N 1 -n 1 -c 2 --cpu_bind=cores ipcontroller --ip='*' --nodb --cluster-id=$cluster_id &
sleep 1
sleep 60
srun -N 238 -n 7600 -c 2 --cpu_bind=cores ipengine --mpi=mpi4py --cluster-id=$cluster_id &
sleep 1
sleep 300
srun -N 1 -n 1 -c 2 --cpu_bind=cores python parallel_optimize.py --config-file-path='data/parallel_optimize_BTSP_CA1_v3_full_config.yaml' --cluster-id=$cluster_id --pop-size=200 --max-iter=50 --path-length=3 --disp --export --output-dir=$SCRATCH/BTSP
