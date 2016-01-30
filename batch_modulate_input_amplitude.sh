#!/bin/bash
#cd $HOME/PycharmProjects/NEURON/
cd $HOME/CA1Sim
qsub -N job2_013016_e2200_i400_subt_modinh0_10 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_subtractive_modinh.py 0 2200 400 0 10'
qsub -N job2_013016_e2200_i400_subt_modinh0_11 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_subtractive_modinh.py 0 2200 400 0 11'
qsub -N job2_013016_e2400_i400_subt_modinh0_12 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_subtractive_modinh.py 0 2400 400 0 12'
qsub -N job2_013016_e2400_i400_subt_modinh0_13 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_subtractive_modinh.py 0 2400 400 0 13'
qsub -N job2_013016_e2600_i400_subt_modinh0_14 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_subtractive_modinh.py 0 2600 400 0 14'
qsub -N job2_013016_e2600_i400_subt_modinh0_15 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_subtractive_modinh.py 0 2600 400 0 15'
qsub -N job2_013016_e1800_i400_subt_modinh0_16 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_subtractive_modinh.py 0 1800 400 0 16'
qsub -N job2_013016_e1800_i400_subt_modinh0_17 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_subtractive_modinh.py 0 1800 400 0 17'
qsub -N job2_013016_e2000_i400_subt_modinh0_18 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_subtractive_modinh.py 0 2000 400 0 18'
qsub -N job2_013016_e2000_i400_subt_modinh0_19 -l haswell=true -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs_subtractive_modinh.py 0 2000 400 0 19'