#!/bin/bash
#cd $HOME/PycharmProjects/NEURON/
cd $HOME/CA1Sim
qsub -N job_012916_e1000_i600_subt_modinh0_0 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs_subtractive_modinh.py 0 1000 600 0 0'
qsub -N job_012916_e1000_i600_subt_modinh0_1 -b y -cwd -V 'python test_poisson_inputs_subtractive_modinh.py 0 1000 600 0 1'
qsub -N job_012916_e1200_i600_subt_modinh0_2 -b y -cwd -V 'python test_poisson_inputs_subtractive_modinh.py 0 1200 600 0 2'
qsub -N job_012916_e1200_i600_subt_modinh0_3 -b y -cwd -V 'python test_poisson_inputs_subtractive_modinh.py 0 1200 600 0 3'
qsub -N job_012916_e1400_i600_subt_modinh0_4 -b y -cwd -V 'python test_poisson_inputs_subtractive_modinh.py 0 1400 600 0 4'
qsub -N job_012916_e1400_i600_subt_modinh0_5 -b y -cwd -V 'python test_poisson_inputs_subtractive_modinh.py 0 1400 600 0 5'
qsub -N job_012916_e1600_i600_subt_modinh0_6 -b y -cwd -V 'python test_poisson_inputs_subtractive_modinh.py 0 1600 600 0 6'
qsub -N job_012916_e1600_i600_subt_modinh0_7 -b y -cwd -V 'python test_poisson_inputs_subtractive_modinh.py 0 1600 600 0 7'
qsub -N job_012916_e1800_i600_subt_modinh0_8 -b y -cwd -V 'python test_poisson_inputs_subtractive_modinh.py 0 1800 600 0 8'
qsub -N job_012916_e1800_i600_subt_modinh0_9 -b y -cwd -V 'python test_poisson_inputs_subtractive_modinh.py 0 1800 600 0 9'