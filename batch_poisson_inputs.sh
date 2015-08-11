#!/bin/bash
#cd $HOME/PycharmProjects/NEURON/
cd $HOME/CA1Sim
qsub -N job_081215_i200_0 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs.py 1 1400 200 0'
qsub -N job_081215_i200_1 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs.py 1 1400 200 1'
qsub -N job_081215_i200_2 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs.py 1 1400 200 2'
qsub -N job_081215_i200_3 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs.py 1 1400 200 3'
qsub -N job_081215_i400_0 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs.py 1 1400 400 0'
qsub -N job_081215_i400_1 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs.py 1 1400 400 1'
qsub -N job_081215_i400_2 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs.py 1 1400 400 2'
qsub -N job_081215_i400_3 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs.py 1 1400 400 3'
qsub -N job_081215_i600_0 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs.py 1 1400 600 0'
qsub -N job_081215_i600_1 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs.py 1 1400 600 1'
qsub -N job_081215_i600_2 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs.py 1 1400 600 2'
qsub -N job_081215_i600_3 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs.py 1 1400 600 3'