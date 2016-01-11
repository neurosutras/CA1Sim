#!/bin/bash
#cd $HOME/PycharmProjects/NEURON/
cd $HOME/CA1Sim
qsub -N job_011116_e1600_i600_no_nmda_inh0_90 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs_no_nmda.py 3 1600 600 0 90'
qsub -N job_011116_e1600_i600_no_nmda_inh0_91 -b y -cwd -V 'python test_poisson_inputs_no_nmda.py 3 1600 600 0 91'
qsub -N job_011116_e1600_i600_no_nmda_inh0_92 -b y -cwd -V 'python test_poisson_inputs_no_nmda.py 3 1600 600 0 92'
qsub -N job_011116_e1600_i600_no_nmda_inh0_93 -b y -cwd -V 'python test_poisson_inputs_no_nmda.py 3 1600 600 0 93'
qsub -N job_011116_e2000_i600_no_nmda_inh0_100 -b y -cwd -V 'python test_poisson_inputs_no_nmda.py 3 2000 600 0 100'
qsub -N job_011116_e2000_i600_no_nmda_inh0_101 -b y -cwd -V 'python test_poisson_inputs_no_nmda.py 3 2000 600 0 101'
qsub -N job_011116_e2000_i600_no_nmda_inh0_102 -b y -cwd -V 'python test_poisson_inputs_no_nmda.py 3 2000 600 0 102'
qsub -N job_011116_e2000_i600_no_nmda_inh0_103 -b y -cwd -V 'python test_poisson_inputs_no_nmda.py 3 2000 600 0 103'
qsub -N job_011116_e2400_i600_no_nmda_inh0_110 -b y -cwd -V 'python test_poisson_inputs_no_nmda.py 3 2400 600 0 110'
qsub -N job_011116_e2400_i600_no_nmda_inh0_111 -b y -cwd -V 'python test_poisson_inputs_no_nmda.py 3 2400 600 0 111'
qsub -N job_011116_e2400_i600_no_nmda_inh0_112 -b y -cwd -V 'python test_poisson_inputs_no_nmda.py 3 2400 600 0 112'
qsub -N job_011116_e2400_i600_no_nmda_inh0_113 -b y -cwd -V 'python test_poisson_inputs_no_nmda.py 3 2400 600 0 113'