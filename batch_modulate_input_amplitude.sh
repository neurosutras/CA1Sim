#!/bin/bash
#cd $HOME/PycharmProjects/NEURON/
cd $HOME/CA1Sim
qsub -N job_013016_e2200_i600_subt_modinh0_0 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_subtractive_modinh.py 0 2200 600 0 0'
qsub -N job_013016_e2200_i600_subt_modinh0_1 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_subtractive_modinh.py 0 2200 600 0 1'
qsub -N job_013016_e2400_i600_subt_modinh0_2 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_subtractive_modinh.py 0 2400 600 0 2'
qsub -N job_013016_e2400_i600_subt_modinh0_3 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_subtractive_modinh.py 0 2400 600 0 3'
qsub -N job_013016_e2600_i600_subt_modinh0_4 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_subtractive_modinh.py 0 2600 600 0 4'
qsub -N job_013016_e2600_i600_subt_modinh0_5 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_subtractive_modinh.py 0 2600 600 0 5'
qsub -N job_013016_e2800_i600_subt_modinh0_6 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_subtractive_modinh.py 0 2800 600 0 6'
qsub -N job_013016_e2800_i600_subt_modinh0_7 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_subtractive_modinh.py 0 2800 600 0 7'
qsub -N job_013016_e3000_i600_subt_modinh0_8 -l haswell=true -b y -cwd -V 'python test_poisson_inputs_subtractive_modinh.py 0 3000 600 0 8'
qsub -N job_013016_e3000_i600_subt_modinh0_9 -l haswell=true -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs_subtractive_modinh.py 0 3000 600 0 9'