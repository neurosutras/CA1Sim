#!/bin/bash
#
#SBATCH -J optimize_leak
#SBATCH -o optimize_leak_062717.%j.o
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
group_sz="$2"
pop_size="$3"
max_iter="$4"
ipcontroller --ip='*' --quiet --cluster-id=$cluster_id &
sleep 60
mpirun -np 300 ipengine --quiet --cluster-id=$cluster_id &
sleep 60
python parallel_optimize_leak.py --cluster-id=$cluster_id --group-size=$group_sz --pop-size=$pop_size --max-iter=$max_iter --optimize --disp