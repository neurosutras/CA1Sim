#!/bin/bash
#SBATCH -J parallel_optimize_leak_comet
#SBATCH -o parallel_optimize_leak_comet_071317.%j.o
#SBATCH -e optimize_leak_comet_071317.%j.e
#SBATCH -N 26 -n 602
#SBATCH -t 00:30:00
##SBATCH --mail-user=graceyng@stanford.edu
#SBATCH --mail-user=aaronmil@stanford.edu
#SBATCH --mail-type=END
#SBATCH --mail-type=BEGIN
#

set -x

cd $HOME/CA1Sim
cluster_id="test_comet"

ibrun -n 1 ipcontroller --ip='*' --cluster-id=$cluster_id &
sleep 1
sleep 60
ibrun -n 600 ipengine --cluster-id=$cluster_id &
sleep 1
sleep 180
ibrun -n 1 python parallel_optimize_leak.py --cluster-id=$cluster_id --group-size=3 --pop-size=200 --max-iter=2 --path-length=2 --disp
