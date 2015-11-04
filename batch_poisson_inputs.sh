#!/bin/bash
#cd $HOME/PycharmProjects/NEURON/
cd $HOME/CA1Sim
qsub -N job_110615_e3600_i400_modamp_0 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs_modulate_amplitude.py 1 3600 400 0'
qsub -N job_110615_e3600_i400_modamp_1 -b y -cwd -V 'python test_poisson_inputs_modulate_amplitude.py 1 3600 400 1'
qsub -N job_110615_e3600_i400_modamp_2 -b y -cwd -V 'python test_poisson_inputs_modulate_amplitude.py 1 3600 400 2'
qsub -N job_110615_e3600_i400_modamp_3 -b y -cwd -V 'python test_poisson_inputs_modulate_amplitude.py 1 3600 400 3'
qsub -N job_110615_e3600_i200_modamp_10 -b y -cwd -V 'python test_poisson_inputs_modulate_amplitude.py 1 3600 200 10'
qsub -N job_110615_e3600_i200_modamp_11 -b y -cwd -V 'python test_poisson_inputs_modulate_amplitude.py 1 3600 200 11'
qsub -N job_110615_e3600_i200_modamp_12 -b y -cwd -V 'python test_poisson_inputs_modulate_amplitude.py 1 3600 200 12'
qsub -N job_110615_e3600_i200_modamp_13 -b y -cwd -V 'python test_poisson_inputs_modulate_amplitude.py 1 3600 200 13'