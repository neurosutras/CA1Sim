#!/bin/bash
#cd $HOME/PycharmProjects/NEURON/
cd $HOME/CA1Sim
qsub -N job_122215_e1600_i600_cal_fb_inh0_0 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs_calibrate_feedback.py 0 1600 600 0 0'
qsub -N job_122215_e1600_i600_cal_fb_inh0_1 -b y -cwd -V 'python test_poisson_inputs_calibrate_feedback.py 0 1600 600 0 1'
qsub -N job_122215_e1600_i600_cal_fb_inh0_2 -b y -cwd -V 'python test_poisson_inputs_calibrate_feedback.py 0 1600 600 0 2'
qsub -N job_122215_e1600_i600_cal_fb_inh0_3 -b y -cwd -V 'python test_poisson_inputs_calibrate_feedback.py 0 1600 600 0 3'
qsub -N job_122215_e1600_i600_cal_fb_inh1_10 -b y -cwd -V 'python test_poisson_inputs_calibrate_feedback.py 0 1600 600 1 10'
qsub -N job_122215_e1600_i600_cal_fb_inh1_11 -b y -cwd -V 'python test_poisson_inputs_calibrate_feedback.py 0 1600 600 1 11'
qsub -N job_122215_e1600_i600_cal_fb_inh1_12 -b y -cwd -V 'python test_poisson_inputs_calibrate_feedback.py 0 1600 600 1 12'
qsub -N job_122215_e1600_i600_cal_fb_inh1_13 -b y -cwd -V 'python test_poisson_inputs_calibrate_feedback.py 0 1600 600 1 13'
qsub -N job_122215_e1600_i600_cal_fb_inh2_20 -b y -cwd -V 'python test_poisson_inputs_calibrate_feedback.py 0 1600 600 2 20'
qsub -N job_122215_e1600_i600_cal_fb_inh2_21 -b y -cwd -V 'python test_poisson_inputs_calibrate_feedback.py 0 1600 600 2 21'
qsub -N job_122215_e1600_i600_cal_fb_inh2_22 -b y -cwd -V 'python test_poisson_inputs_calibrate_feedback.py 0 1600 600 2 22'
qsub -N job_122215_e1600_i600_cal_fb_inh2_23 -b y -cwd -V 'python test_poisson_inputs_calibrate_feedback.py 0 1600 600 2 23'