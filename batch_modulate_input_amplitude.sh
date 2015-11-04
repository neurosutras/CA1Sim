#!/bin/bash
#cd $HOME/PycharmProjects/NEURON/
cd $HOME/CA1Sim
qsub -N job_110415_e6000_i400_mod_amp_10 -l haswell=true -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs_modulate_amplitude.py 1 6000 400 10'
qsub -N job_110415_e6000_i400_mod_amp_11 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_modulate_amplitude.py 1 6000 400 11'
qsub -N job_110415_e6000_i400_mod_amp_12 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_modulate_amplitude.py 1 6000 400 12'
qsub -N job_110415_e6000_i400_mod_amp_13 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_modulate_amplitude.py 1 6000 400 13'
qsub -N job_110415_e6000_i400_mod_num_20 -l haswell=true -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs_modulate_number.py 1 6000 400 20'
qsub -N job_110415_e6000_i400_mod_num_21 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_modulate_number.py 1 6000 400 21'
qsub -N job_110415_e6000_i400_mod_num_22 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_modulate_number.py 1 6000 400 22'
qsub -N job_110415_e6000_i400_mod_num_23 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_modulate_number.py 1 6000 400 23'