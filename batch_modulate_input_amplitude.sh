#!/bin/bash
#cd $HOME/PycharmProjects/NEURON/
cd $HOME/CA1Sim
qsub -N job_121315_e1550_i600_cal_ec_inh0_7 -l haswell=true -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs_calibrate_ECIII.py 0 1550 600 0 7'
qsub -N job_121315_e1550_i600_cal_ec_inh0_8 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_calibrate_ECIII.py 0 1550 600 0 8'
qsub -N job_121315_e1550_i600_cal_ec_inh0_9 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_calibrate_ECIII.py 0 1550 600 0 9'
qsub -N job_121315_e1550_i600_cal_ec_inh0_10 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_calibrate_ECIII.py 0 1550 600 0 10'
qsub -N job_121315_e1550_i600_cal_ec_inh2_11 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_calibrate_ECIII.py 0 1550 600 2 11'
qsub -N job_121315_e1550_i600_cal_ec_inh2_12 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_calibrate_ECIII.py 0 1550 600 2 12'
qsub -N job_121315_e1550_i600_cal_ec_inh2_13 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_calibrate_ECIII.py 0 1550 600 2 13'
qsub -N job_121315_e1550_i600_cal_ec_inh2_14 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_calibrate_ECIII.py 0 1550 600 2 14'