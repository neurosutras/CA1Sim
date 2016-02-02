#!/bin/bash
#cd $HOME/PycharmProjects/NEURON/
cd $HOME/CA1Sim
qsub -N job2_020216_e3400_i400_cal_pr1_0 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_calibrate_pr.py 1 3400 400 0 0'
qsub -N job2_020216_e3400_i400_cal_pr1_1 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_calibrate_pr.py 1 3400 400 0 1'
qsub -N job2_020216_e3600_i400_cal_pr1_2 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_calibrate_pr.py 1 3600 400 0 2'
qsub -N job2_020216_e3600_i400_cal_pr1_3 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_calibrate_pr.py 1 3600 400 0 3'
qsub -N job2_020216_e3800_i500_cal_pr1_10 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_calibrate_pr.py 1 3800 500 0 10'
qsub -N job2_020216_e3800_i500_cal_pr1_11 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_calibrate_pr.py 1 3800 500 0 11'
qsub -N job2_020216_e4000_i500_cal_pr1_12 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_calibrate_pr.py 1 4000 500 0 12'
qsub -N job2_020216_e4000_i500_cal_pr1_13 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_calibrate_pr.py 1 4000 500 0 13'
qsub -N job2_020216_e4400_i600_cal_pr1_20 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_calibrate_pr.py 1 4400 600 0 20'
qsub -N job2_020216_e4400_i600_cal_pr1_21 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_calibrate_pr.py 1 4400 600 0 21'
qsub -N job2_020216_e4600_i600_cal_pr1_22 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_calibrate_pr.py 1 4600 600 0 22'
qsub -N job2_020216_e4600_i600_cal_pr1_23 -l haswell=true -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs_calibrate_pr.py 1 4600 600 0 23'