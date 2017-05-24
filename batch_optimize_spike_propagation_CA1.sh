#!/bin/bash
#
#SBATCH -J optimize_spike_propagation_CA1
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH --mem-per-cpu=4G
##SBATCH --mail-user=graceyng@stanford.edu
#SBATCH --mail-user=aaronmil@stanford.edu
#SBATCH --mail-type=END
#SBATCH --mail-type=BEGIN
#
set -x

cd $HOME/CA1Sim
python optimize_spike_propagation_CA1Pyr.py
