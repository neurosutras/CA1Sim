#!/bin/bash
#cd $HOME/PycharmProjects/NEURON/
cd $HOME/CA1Sim
qsub -N job_121615_e1700_i600_cal_fb_inh0_24 -l haswell=true -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs_calibrate_feedback.py 0 1700 600 0 24'
qsub -N job_121615_e1700_i600_cal_fb_inh0_25 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_calibrate_feedback.py 0 1700 600 0 25'
qsub -N job_121615_e1700_i600_cal_fb_inh0_26 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_calibrate_feedback.py 0 1700 600 0 26'
qsub -N job_121615_e1700_i600_cal_fb_inh0_27 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_calibrate_feedback.py 0 1700 600 0 27'
qsub -N job_121615_e1700_i600_cal_fb_inh1_44 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_calibrate_feedback.py 0 1700 600 1 44'
qsub -N job_121615_e1700_i600_cal_fb_inh1_45 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_calibrate_feedback.py 0 1700 600 1 45'
qsub -N job_121615_e1700_i600_cal_fb_inh1_46 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_calibrate_feedback.py 0 1700 600 1 46'
qsub -N job_121615_e1700_i600_cal_fb_inh1_47 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_calibrate_feedback.py 0 1700 600 1 47'
qsub -N job_121615_e1700_i600_cal_fb_inh2_34 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_calibrate_feedback.py 0 1700 600 2 34'
qsub -N job_121615_e1700_i600_cal_fb_inh2_35 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_calibrate_feedback.py 0 1700 600 2 35'
qsub -N job_121615_e1700_i600_cal_fb_inh2_36 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_calibrate_feedback.py 0 1700 600 2 36'
qsub -N job_121615_e1700_i600_cal_fb_inh2_37 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_calibrate_feedback.py 0 1700 600 2 37'