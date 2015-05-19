#!/bin/bash
#cd $HOME/.ipython/
#cd $HOME/PycharmProjects/NEURON/
cd $HOME/CA1Sim
#ipcluster start -n 8 --ip='*' --profile=mpi &
ipcluster start -n 8 --ip='*' --profile-dir=$HOME/.ipython/profile_default --engines=MPIEngineSetLauncher &
#ipcluster start -n 32 --ip='*' --profile=mpi &
sleep 60
#cd $HOME/CA1Sim/
python testMPI_controller.py
#cd $HOME/.ipython/
#ipcluster stop --profile=mpi
ipcluster stop
