#!/bin/bash
#PBS -l nodes=39:ppn=16:xe
#PBS -q debug
#PBS -l walltime=00:30:00
#PBS -A bafv
#PBS -N parallel_optimize_leak_bw_071317
#PBS -e parallel_optimize_leak_bw_071317.$PBS_JOBID.err
#PBS -o parallel_optimize_leak_bw_071317.$PBS_JOBID.out
#PBS -m bea
#PBS -M aaronmil@stanford.edu
#PBS -W umask=0027

module swap PrgEnv-cray PrgEnv-gnu
module load cray-hdf5-parallel
module load bwpy
module load bwpy-mpi

export PI_HOME=/projects/sciteam/baef
# export LD_LIBRARY_PATH=/sw/bw/bwpy/0.3.0/python-mpi/usr/lib:/sw/bw/bwpy/0.3.0/python-single/usr/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$PI_HOME/bin/nrn/lib/python:$PI_HOME/site-packages:$PYTHONPATH
export PATH=$PI_HOME/bin/nrn/x86_64/bin:$PATH
cd $HOME/CA1Sim

cluster_id="leak"
#profile="bw"
profile='mpi'

set -x

# max number of engines = (nodes - 1) * ppn - 1
ipcluster start -n 607 --daemon --ip='*' --cluster-id=$cluster_id --profile=$profile &
sleep 1
sleep 1
aprun -n 1 python parallel_optimize_leak.py --cluster-id=$cluster_id --profile=$profile --group-size=3 --pop-size=200 --path-length=2 --max-iter=2 --disp --sleep
ipcluster stop --cluster-id=$cluster_id --profile=$profile
