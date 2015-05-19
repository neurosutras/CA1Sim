#!/bin/bash
cd $HOME/.ipython/
ipcluster start -n 8 --ip='*' --profile=mpi &
#ipcluster start -n 32 --ip='*' --profile=mpi &
#PID=$!
sleep 60
#cd $HOME/PycharmProjects/NEURON/
cd $HOME/CA1Sim
python parallel_branch_cooperativity_controller.py
#sleep 10
#kill -INT $PID
cd $HOME/.ipython/
ipcluster stop --profile=mpi
