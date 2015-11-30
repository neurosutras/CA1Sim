#!/bin/bash
#cd $HOME/PycharmProjects/NEURON/
cd $HOME/CA1Sim
qsub -N job_113015_e1250_i300_mod_amp_inh0_20 -l haswell=true -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs_modulate_amplitude.py 2 1250 300 0 20'
qsub -N job_113015_e1250_i300_mod_amp_inh0_21 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_modulate_amplitude.py 1 1250 300 0 21'
qsub -N job_113015_e1250_i300_mod_amp_inh0_22 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_modulate_amplitude.py 1 1250 300 0 22'
qsub -N job_113015_e1250_i300_mod_amp_inh0_23 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_modulate_amplitude.py 1 1250 300 0 23'
qsub -N job_113015_e1250_i300_mod_amp_inh1_30 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_modulate_amplitude.py 1 1250 300 1 30'
qsub -N job_113015_e1250_i300_mod_amp_inh1_31 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_modulate_amplitude.py 1 1250 300 1 31'
qsub -N job_113015_e1250_i300_mod_amp_inh1_32 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_modulate_amplitude.py 1 1250 300 1 32'
qsub -N job_113015_e1250_i300_mod_amp_inh1_33 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_modulate_amplitude.py 1 1250 300 1 33'
qsub -N job_113015_e1250_i300_mod_amp_inh2_40 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_modulate_amplitude.py 1 1250 300 2 40'
qsub -N job_113015_e1250_i300_mod_amp_inh2_41 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_modulate_amplitude.py 1 1250 300 2 41'
qsub -N job_113015_e1250_i300_mod_amp_inh2_42 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_modulate_amplitude.py 1 1250 300 2 42'
qsub -N job_113015_e1250_i300_mod_amp_inh2_43 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_modulate_amplitude.py 1 1250 300 2 43'