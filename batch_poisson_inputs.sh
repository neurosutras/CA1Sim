#!/bin/bash
#cd $HOME/PycharmProjects/NEURON/
cd $HOME/CA1Sim
qsub -N job_100115_e1400_i0_0 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs.py 1 1400 0 0'
qsub -N job_100115_e1400_i0_1 -b y -cwd -V 'python test_poisson_inputs.py 1 1400 0 1'
qsub -N job_100115_e1400_i0_2 -b y -cwd -V 'python test_poisson_inputs.py 1 1400 0 2'
qsub -N job_100115_e1400_i0_3 -b y -cwd -V 'python test_poisson_inputs.py 1 1400 0 3'
qsub -N job_100115_e1600_i0_0 -b y -cwd -V 'python test_poisson_inputs.py 1 1600 0 0'
qsub -N job_100115_e1600_i0_1 -b y -cwd -V 'python test_poisson_inputs.py 1 1600 0 1'
qsub -N job_100115_e1600_i0_2 -b y -cwd -V 'python test_poisson_inputs.py 1 1600 0 2'
qsub -N job_100115_e1600_i0_3 -b y -cwd -V 'python test_poisson_inputs.py 1 1600 0 3'
qsub -N job_100115_e1800_i0_0 -b y -cwd -V 'python test_poisson_inputs.py 1 1800 0 0'
qsub -N job_100115_e1800_i0_1 -b y -cwd -V 'python test_poisson_inputs.py 1 1800 0 1'
qsub -N job_100115_e1800_i0_2 -b y -cwd -V 'python test_poisson_inputs.py 1 1800 0 2'
qsub -N job_100115_e1800_i0_3 -b y -cwd -V 'python test_poisson_inputs.py 1 1800 0 3'
qsub -N job_100115_e2000_i0_0 -b y -cwd -V 'python test_poisson_inputs.py 1 2000 0 0'
qsub -N job_100115_e2000_i0_1 -b y -cwd -V 'python test_poisson_inputs.py 1 2000 0 1'
qsub -N job_100115_e2000_i0_2 -b y -cwd -V 'python test_poisson_inputs.py 1 2000 0 2'
qsub -N job_100115_e2000_i0_3 -b y -cwd -V 'python test_poisson_inputs.py 1 2000 0 3'