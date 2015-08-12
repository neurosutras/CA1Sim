#!/bin/bash
#cd $HOME/PycharmProjects/NEURON/
cd $HOME/CA1Sim
qsub -N job_081215_i400_mod_amp_0 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs_modulate_amplitude.py 1 1400 400 0'
qsub -N job_081215_i400_mod_amp_1 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs_modulate_amplitude.py 1 1400 400 1'
qsub -N job_081215_i400_mod_amp_2 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs_modulate_amplitude.py 1 1400 400 2'
qsub -N job_081215_i400_mod_amp_3 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs_modulate_amplitude.py 1 1400 400 3'