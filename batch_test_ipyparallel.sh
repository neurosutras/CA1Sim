#!/bin/bash
#cd $HOME/PycharmProjects/NEURON/
#SBATCH -N 2 -n 31 -J 110216_test_ipyparallel -o 110216.o --mem-per-cpu 1G
cd $HOME/CA1Sim
cluster_id="$1"
ipcluster start -n 31 --profile-dir=$HOME/.ipython/profile_default --cluster-id=$cluster_id &
sleep 180
ipython parallel_test_ipyparallel_controller.py $cluster_id
ipcluster stop --profile-dir=$HOME/.ipython/profile_default --cluster-id=$cluster_id