#!/bin/bash
#cd $HOME/PycharmProjects/NEURON/
cd $HOME/CA1Sim
qsub -N job_011116_e2000_i600_no_nmda_inh0_101 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs_no_nmda.py 3 2000 600 0 101'
qsub -N job_011116_e2000_i600_no_nmda_inh0_102 -b y -cwd -V 'python test_poisson_inputs_no_nmda.py 3 2000 600 0 102'
qsub -N job_011116_e2000_i600_no_nmda_inh0_103 -b y -cwd -V 'python test_poisson_inputs_no_nmda.py 3 2000 600 0 103'
qsub -N job_011116_e2400_i600_no_nmda_inh0_110 -b y -cwd -V 'python test_poisson_inputs_no_nmda.py 3 2400 600 0 110'
qsub -N job_011116_e2400_i600_no_nmda_inh0_111 -b y -cwd -V 'python test_poisson_inputs_no_nmda.py 3 2400 600 0 111'