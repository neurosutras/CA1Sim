#!/bin/bash
cd $HOME/CA1Sim
ipcluster start -n 32 --engines=MPI &
sleep 60
python parallel_branch_cooperativity_controller.py
sleep 10
ipcluster stop
wait
