#!/bin/bash
#cd $HOME/PycharmProjects/NEURON/
cd $HOME/CA1Sim
qsub -N job_011816_e3600_i600_no_nmda_inh1_140 -l haswell=true -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs_no_nmda.py 3 3600 600 1 140'
qsub -N job_011816_e3600_i600_no_nmda_inh1_141 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_no_nmda.py 3 3600 600 1 141'
qsub -N job_011816_e3600_i600_no_nmda_inh1_142 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_no_nmda.py 3 3600 600 1 142'
qsub -N job_011816_e3600_i600_no_nmda_inh1_143 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_no_nmda.py 3 3600 600 1 143'
qsub -N job_011816_e3600_i600_no_nmda_inh2_150 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_no_nmda.py 3 3600 600 2 150'
qsub -N job_011816_e3600_i600_no_nmda_inh2_151 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_no_nmda.py 3 3600 600 2 151'
qsub -N job_011816_e3600_i600_no_nmda_inh2_152 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_no_nmda.py 3 3600 600 2 152'
qsub -N job_011816_e3600_i600_no_nmda_inh2_153 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_no_nmda.py 3 3600 600 2 153'