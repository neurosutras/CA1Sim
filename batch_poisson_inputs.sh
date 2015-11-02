#!/bin/bash
#cd $HOME/PycharmProjects/NEURON/
cd $HOME/CA1Sim
qsub -N job_110215_e6000_i0_0 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs.py 1 6000 0 0'
qsub -N job_110215_e6000_i0_1 -b y -cwd -V 'python test_poisson_inputs.py 1 6000 0 1'
qsub -N job_110215_e6000_i0_2 -b y -cwd -V 'python test_poisson_inputs.py 1 6000 0 2'
qsub -N job_110215_e6000_i0_3 -b y -cwd -V 'python test_poisson_inputs.py 1 6000 0 3'
qsub -N job_110215_e4800_i0_0 -b y -cwd -V 'python test_poisson_inputs.py 1 4800 0 0'
qsub -N job_110215_e4800_i0_1 -b y -cwd -V 'python test_poisson_inputs.py 1 4800 0 1'
qsub -N job_110215_e4800_i0_2 -b y -cwd -V 'python test_poisson_inputs.py 1 4800 0 2'
qsub -N job_110215_e4800_i0_3 -b y -cwd -V 'python test_poisson_inputs.py 1 4800 0 3'
qsub -N job_110215_e3600_i0_0 -b y -cwd -V 'python test_poisson_inputs.py 1 3600 0 0'
qsub -N job_110215_e3600_i0_1 -b y -cwd -V 'python test_poisson_inputs.py 1 3600 0 1'
qsub -N job_110215_e3600_i0_2 -b y -cwd -V 'python test_poisson_inputs.py 1 3600 0 2'
qsub -N job_110215_e3600_i0_3 -b y -cwd -V 'python test_poisson_inputs.py 1 3600 0 3'