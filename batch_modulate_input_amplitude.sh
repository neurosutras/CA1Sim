#!/bin/bash
#cd $HOME/PycharmProjects/NEURON/
cd $HOME/CA1Sim
qsub -N job_121715_e1700_i550_cal_fb_inh0_50 -l haswell=true -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs_calibrate_feedback.py 0 1700 550 0 50'
qsub -N job_121715_e1700_i550_cal_fb_inh0_51 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_calibrate_feedback.py 0 1700 550 0 51'
qsub -N job_121715_e1700_i550_cal_fb_inh0_52 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_calibrate_feedback.py 0 1700 550 0 52'
qsub -N job_121715_e1700_i550_cal_fb_inh0_53 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_calibrate_feedback.py 0 1700 550 0 53'
qsub -N job_121715_e1700_i550_cal_fb_inh1_60 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_calibrate_feedback.py 0 1700 550 1 60'
qsub -N job_121715_e1700_i550_cal_fb_inh1_61 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_calibrate_feedback.py 0 1700 550 1 61'
qsub -N job_121715_e1700_i550_cal_fb_inh1_62 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_calibrate_feedback.py 0 1700 550 1 62'
qsub -N job_121715_e1700_i550_cal_fb_inh1_63 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_calibrate_feedback.py 0 1700 550 1 63'
qsub -N job_121715_e1700_i550_cal_fb_inh2_70 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_calibrate_feedback.py 0 1700 550 2 70'
qsub -N job_121715_e1700_i550_cal_fb_inh2_71 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_calibrate_feedback.py 0 1700 550 2 71'
qsub -N job_121715_e1700_i550_cal_fb_inh2_72 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_calibrate_feedback.py 0 1700 550 2 72'
qsub -N job_121715_e1700_i550_cal_fb_inh2_73 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_calibrate_feedback.py 0 1700 550 2 73'