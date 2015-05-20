#!/bin/bash
#cd $HOME/.ipython/
cd $HOME/PycharmProjects/NEURON/
#cd $HOME/CA1Sim
#ipcluster start -n 8 --ip='*' --profile=mpi &
#ipcluster mpiexec --n=8 --ip='*' --mpi=mpi4py --profile-dir=$HOME/.ipython/profile_default --engines=MPIEngineSetLauncher &
#ipcluster start -n 32 --ip='*' --profile=mpi &
ipcontroller --profile-dir=$HOME/.ipython/profile_default &
sleep 30
mpirun -n 8 ipengine --profile-dir=$HOME/.ipython/profile_default &
sleep 60
#cd $HOME/CA1Sim/
python testMPI_controller.py
#cd $HOME/.ipython/
#ipcluster stop --profile=mpi
ipcluster stop --profile-dir=$HOME/.ipython/profile_default
