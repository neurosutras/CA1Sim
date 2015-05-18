#!/bin/bash
ipcluster start -n 32 --engines=MPI &
PID=$!
sleep 60
python parallel_branch_cooperativity_controller.py
sleep 10
ipcluster stop
wait
