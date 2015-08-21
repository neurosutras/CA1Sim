#!/bin/bash
#cd $HOME/PycharmProjects/NEURON/
cd $HOME/CA1Sim
qsub -N job_082115_i400_mod_num_0 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs_modulate_number.py 1 1400 400 0'
qsub -N job_082115_i400_mod_num_1 -b y -cwd -V 'python test_poisson_inputs_modulate_number.py 1 1400 400 1'
qsub -N job_082115_i400_mod_num_2 -b y -cwd -V 'python test_poisson_inputs_modulate_number.py 1 1400 400 2'
qsub -N job_082115_i400_mod_num_3 -b y -cwd -V 'python test_poisson_inputs_modulate_number.py 1 1400 400 3'
qsub -N job_082115_i400_mod_num_no_nmda_0 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs_modulate_number_no_nmda.py 1 1400 400 0'
qsub -N job_082115_i400_mod_num_no_nmda_1 -b y -cwd -V 'python test_poisson_inputs_modulate_number_no_nmda.py 1 1400 400 1'
qsub -N job_082115_i400_mod_num_no_nmda_2 -b y -cwd -V 'python test_poisson_inputs_modulate_number_no_nmda.py 1 1400 400 2'
qsub -N job_082115_i400_mod_num_no_nmda_3 -b y -cwd -V 'python test_poisson_inputs_modulate_number_no_nmda.py 1 1400 400 3'