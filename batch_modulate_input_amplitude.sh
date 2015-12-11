#!/bin/bash
#cd $HOME/PycharmProjects/NEURON/
cd $HOME/CA1Sim
qsub -N job_121115_e1550_i550_cal_ec_inh0_3 -l haswell=true -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs_calibrate_ECIII.py 0 1550 550 0 3'
qsub -N job_121115_e1550_i550_cal_ec_inh0_4 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_calibrate_ECIII.py 0 1550 550 0 4'
qsub -N job_121115_e1550_i550_cal_ec_inh0_5 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_calibrate_ECIII.py 0 1550 550 0 5'
qsub -N job_121115_e1550_i550_cal_ec_inh0_6 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_calibrate_ECIII.py 0 1550 550 0 6'