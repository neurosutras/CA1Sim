#!/bin/bash
#cd $HOME/PycharmProjects/NEURON/
cd $HOME/CA1Sim
qsub -N job_081015_trial0 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs.py 1 1400 0'
qsub -N job_081015_trial1 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs.py 1 1400 1'
qsub -N job_081015_trial2 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs.py 1 1400 2'
qsub -N job_081015_trial3 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs.py 1 1400 3'