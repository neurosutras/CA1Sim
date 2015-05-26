#!/bin/bash
#cd $HOME/PycharmProjects/NEURON/
cd $HOME/CA1Sim
ipcluster start -n 25 --profile-dir=$HOME/.ipython/profile_default &
#ipcluster start -n 8 --profile-dir=$HOME/.ipython/profile_default &
#ipcluster start -n 16 --profile-dir=$HOME/.ipython/profile_default &
sleep 60
python parallel_optimize_branch_cooperativity_nmda_controller.py
ipcluster stop --profile-dir=$HOME/.ipython/profile_default
