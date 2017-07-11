#!/bin/bash
#
#SBATCH -J optimize_leak
#SBATCH -o optimize_leak_071017.%j.o
#SBATCH -e optimize_leak_071017.%j.e
#SBATCH -N 3 -n 72
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
ibrun -n 70 ipengine --cluster-id=$cluster_id &
sleep 1
sleep 120
ibrun -n 1 python parallel_optimize_leak.py --cluster-id=$cluster_id --group-size=3 --pop-size=23 --max-iter=2 --optimize --disp
