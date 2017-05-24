#!/bin/bash
#
#SBATCH -J optimize_BTSP_kinetic_rule
#SBATCH -o optimize_BTSP_kinetic_rule_052017.%j.o
#SBATCH -n 31
#SBATCH -t 48:00:00
#SBATCH --mem-per-cpu=4G
##SBATCH --mail-user=graceyng@stanford.edu
#SBATCH --mail-user=aaronmil@stanford.edu
#SBATCH --mail-type=END
#SBATCH --mail-type=BEGIN
#
set -x

cd $HOME/CA1Sim
cluster_id="$1"
null_cell_id=1
null_label=cell
ipcontroller --ip='*' --quiet --cluster-id=$cluster_id &
sleep 60
mpirun -np 30 ipengine --quiet --cluster-id=$cluster_id &
sleep 60
python parallel_optimize_BTSP_kinetic_rule_controller_052017.py $null_cell_id $null_label $cluster_id
