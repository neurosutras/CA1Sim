#!/bin/bash
#cd $HOME/PycharmProjects/NEURON/
cd $HOME/CA1Sim
qsub -N job_080715_trial2 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs.py 1 1400 2'
qsub -N job_080715_trial3 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs.py 1 1400 3'
qsub -N job_080715_trial4 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs.py 1 1400 4'
qsub -N job_080715_trial5 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs.py 1 1400 5'
qsub -N job_080715_trial6 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs.py 1 1400 6'
qsub -N job_080715_trial7 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs.py 1 1400 7'
qsub -N job_080715_trial8 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs.py 1 1400 8'
qsub -N job_080715_trial9 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs.py 1 1400 9'