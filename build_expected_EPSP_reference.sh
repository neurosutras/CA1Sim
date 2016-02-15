#!/bin/bash
cd $HOME/CA1Sim
cluster_id="$1"
seed="$2"
ipcluster start -n 15 --profile-dir=$HOME/.ipython/profile_default --cluster-id=$cluster_id &
sleep 180
python build_expected_EPSP_reference_controller.py $cluster_id $seed
ipcluster stop --profile-dir=$HOME/.ipython/profile_default --cluster-id=$cluster_id