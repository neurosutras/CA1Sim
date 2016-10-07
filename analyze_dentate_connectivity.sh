#!/bin/bash
#cd $HOME/PycharmProjects/NEURON/
cd $HOME/CA1Sim_dev
cluster_id="$1"
ipcluster start -n 15 --profile-dir=$HOME/.ipython/profile_default --cluster-id=$cluster_id &
sleep 180
ipython analyze_dentate_connectivity_controller.py $cluster_id
ipcluster stop --profile-dir=$HOME/.ipython/profile_default --cluster-id=$cluster_id