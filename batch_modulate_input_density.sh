#!/bin/bash
cd $HOME/CA1Sim
qsub -N job_022016_e2900_i500_density_0_1.4_0 -b y -cwd -V 'python test_poisson_inputs_modulate_density.py 20 2900 500 0 1.4 0'
qsub -N job_022016_e2900_i500_density_0_1.4_1 -b y -cwd -V 'python test_poisson_inputs_modulate_density.py 20 2900 500 0 1.4 1'
qsub -N job_022016_e2900_i500_density_0_1.5_2 -b y -cwd -V 'python test_poisson_inputs_modulate_density.py 20 2900 500 0 1.5 2'
qsub -N job_022016_e2900_i500_density_0_1.5_3 -b y -cwd -V 'python test_poisson_inputs_modulate_density.py 20 2900 500 0 1.5 3'
qsub -N job_022016_e2900_i500_density_0_1.6_4 -b y -cwd -V 'python test_poisson_inputs_modulate_density.py 20 2900 500 0 1.6 4'
qsub -N job_022016_e2900_i500_density_0_1.6_5 -b y -cwd -V 'python test_poisson_inputs_modulate_density.py 20 2900 500 0 1.6 5'
qsub -N job_022016_e2900_i500_density_0_1.7_6 -b y -cwd -V 'python test_poisson_inputs_modulate_density.py 20 2900 500 0 1.7 6'
qsub -N job_022016_e2900_i500_density_0_1.7_7 -b y -cwd -V 'python test_poisson_inputs_modulate_density.py 20 2900 500 0 1.7 7'
qsub -N job_022016_e2900_i500_density_0_1.8_8 -b y -cwd -V 'python test_poisson_inputs_modulate_density.py 20 2900 500 0 1.8 8'
qsub -N job_022016_e2900_i500_density_0_1.8_9 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs_modulate_density.py 20 2900 500 0 1.8 9'
