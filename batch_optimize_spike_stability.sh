#!/bin/bash
#
#SBATCH -J optimize_spike_stability
#SBATCH -o optimize_spike_stability_041917.%j.o
#SBATCH -n 3
#SBATCH -t 36:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --mail-user=graceyng@stanford.edu
#SBATCH --mail-type=END
#SBATCH --mail-type=BEGIN
#
set -x

cd $HOME/CA1Sim
spines="$1"
mech_filename="$2"
cluster_id="$3"
ipcontroller --ip='*' --quiet --cluster-id=$cluster_id &
sleep 60
mpirun -np 3 ipengine --quiet --cluster-id=$cluster_id &
sleep 60
python parallel_optimize_spike_stability_controller.py $spines "$mech_filename" $cluster_id
