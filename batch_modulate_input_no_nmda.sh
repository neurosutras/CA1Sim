#!/bin/bash
#cd $HOME/PycharmProjects/NEURON/
cd $HOME/CA1Sim
qsub -N job_011516_e3200_i600_no_nmda_inh0_120 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs_no_nmda.py 3 3200 600 0 120'
qsub -N job_011516_e3200_i600_no_nmda_inh0_121 -b y -cwd -V 'python test_poisson_inputs_no_nmda.py 3 3200 600 0 121'
qsub -N job_011516_e3200_i600_no_nmda_inh0_122 -b y -cwd -V 'python test_poisson_inputs_no_nmda.py 3 3200 600 0 122'
qsub -N job_011516_e3200_i600_no_nmda_inh0_123 -b y -cwd -V 'python test_poisson_inputs_no_nmda.py 3 3200 600 0 123'
qsub -N job_011516_e3600_i600_no_nmda_inh0_130 -b y -cwd -V 'python test_poisson_inputs_no_nmda.py 3 3600 600 0 130'
qsub -N job_011516_e3600_i600_no_nmda_inh0_131 -b y -cwd -V 'python test_poisson_inputs_no_nmda.py 3 3600 600 0 131'
qsub -N job_011516_e3600_i600_no_nmda_inh0_132 -b y -cwd -V 'python test_poisson_inputs_no_nmda.py 3 3600 600 0 132'
qsub -N job_011516_e3600_i600_no_nmda_inh0_133 -b y -cwd -V 'python test_poisson_inputs_no_nmda.py 3 3600 600 0 133'