#!/bin/bash
#cd $HOME/PycharmProjects/NEURON/
cd $HOME/CA1Sim
qsub -N job_020216_e1700_i400_subt_modinh0_0 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_subtractive_modinh.py 0 1700 400 0 0'
qsub -N job_020216_e1700_i400_subt_modinh0_1 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_subtractive_modinh.py 0 1700 400 0 1'
qsub -N job_020216_e1700_i400_subt_modinh0_2 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_subtractive_modinh.py 0 1700 400 0 2'
qsub -N job_020216_e1700_i400_subt_modinh0_3 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_subtractive_modinh.py 0 1700 400 0 3'
qsub -N job_020216_e1700_i400_subt_modinh0_6 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_subtractive_modinh.py 0 1700 400 0 6'
qsub -N job_020216_e1700_i400_subt_modinh0_7 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_subtractive_modinh.py 0 1700 400 0 7'
qsub -N job_020216_e1700_i400_subt_modinh0_8 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_subtractive_modinh.py 0 1700 400 0 8'
qsub -N job_020216_e1700_i400_subt_modinh0_9 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_subtractive_modinh.py 0 1700 400 0 9'
qsub -N job_020216_e1900_i500_subt_modinh0_10 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_subtractive_modinh.py 0 1900 500 0 10'
qsub -N job_020216_e1900_i500_subt_modinh0_11 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_subtractive_modinh.py 0 1900 500 0 11'
qsub -N job_020216_e1900_i500_subt_modinh0_14 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_subtractive_modinh.py 0 1900 500 0 14'
qsub -N job_020216_e1900_i500_subt_modinh0_15 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_subtractive_modinh.py 0 1900 500 0 15'
qsub -N job_020216_e1900_i500_subt_modinh0_16 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_subtractive_modinh.py 0 1900 500 0 16'
qsub -N job_020216_e1900_i500_subt_modinh0_17 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_subtractive_modinh.py 0 1900 500 0 17'
qsub -N job_020216_e1900_i500_subt_modinh0_18 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_subtractive_modinh.py 0 1900 500 0 18'
qsub -N job_020216_e1900_i500_subt_modinh0_19 -l haswell=true -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs_subtractive_modinh.py 0 1900 500 0 19'