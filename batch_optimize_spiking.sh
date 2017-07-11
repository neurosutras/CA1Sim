#!/bin/bash
#
#SBATCH -J optimize_spiking
#SBATCH -o optimize_spiking_071117.%j.o
#SBATCH -n 302
#SBATCH -t 36:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --mail-user=graceyng@stanford.edu
##SBATCH --mail-user=aaronmil@stanford.edu
#SBATCH --mail-type=END
#SBATCH --mail-type=BEGIN
#
set -x

cd $HOME/CA1Sim
cluster_id="$1"
group_sz1="$2"
group_sz2="$3"
pop_size="$4"
max_gen="$5"
ipcontroller --ip='*' --quiet --cluster-id=$cluster_id &
sleep 60
mpirun -np 300 ipengine --quiet --cluster-id=$cluster_id &
sleep 60
python parallel_optimize_spiking.py --cluster-id=$cluster_id --group-sizes $group_sz1 $group_sz2 --pop-size=$pop_size --max-gens=$max_gen --optimize --disp