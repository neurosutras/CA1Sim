#!/bin/bash
#cd $HOME/PycharmProjects/NEURON/
cd $HOME/CA1Sim
cluster_id="$1"
ipcluster start -n 31 --profile-dir=$HOME/.ipython/profile_default --cluster-id=$cluster_id &
sleep 180
python parallel_clustered_branch_cooperativity_no_nmda_controller.py $cluster_id
ipcluster stop --profile-dir=$HOME/.ipython/profile_default --cluster-id=$cluster_id