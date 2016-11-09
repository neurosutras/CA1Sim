#!/bin/bash
#cd $HOME/PycharmProjects/NEURON/
#SBATCH -N 1 -n 31 -J 110116_rinp -o 110116.o --mem 8G
cd $HOME/CA1Sim
spines="$1"
cluster_id="$2"
ipcluster start -n 31 --profile-dir=$HOME/.ipython/profile_default --cluster-id=$cluster_id &
sleep 180
ipython parallel_rinp_controller.py $spines $cluster_id
ipcluster stop --profile-dir=$HOME/.ipython/profile_default --cluster-id=$cluster_id