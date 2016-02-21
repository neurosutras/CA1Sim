#!/bin/bash
cd $HOME/CA1Sim
qsub -N job_022016_e2900_i500_density_0_0.4_0 -b y -cwd -V 'python test_poisson_inputs_modulate_density.py 20 2900 500 0 0.4 0'
qsub -N job_022016_e2900_i500_density_0_0.4_1 -b y -cwd -V 'python test_poisson_inputs_modulate_density.py 20 2900 500 0 0.4 1'
qsub -N job_022016_e2900_i500_density_0_0.5_2 -b y -cwd -V 'python test_poisson_inputs_modulate_density.py 20 2900 500 0 0.5 2'
qsub -N job_022016_e2900_i500_density_0_0.5_3 -b y -cwd -V 'python test_poisson_inputs_modulate_density.py 20 2900 500 0 0.5 3'
qsub -N job_022016_e2900_i500_density_0_0.6_4 -b y -cwd -V 'python test_poisson_inputs_modulate_density.py 20 2900 500 0 0.6 4'
qsub -N job_022016_e2900_i500_density_0_0.6_5 -b y -cwd -V 'python test_poisson_inputs_modulate_density.py 20 2900 500 0 0.6 5'
qsub -N job_022016_e2900_i500_density_0_0.7_6 -b y -cwd -V 'python test_poisson_inputs_modulate_density.py 20 2900 500 0 0.7 6'
qsub -N job_022016_e2900_i500_density_0_0.7_7 -b y -cwd -V 'python test_poisson_inputs_modulate_density.py 20 2900 500 0 0.7 7'
qsub -N job_022016_e2900_i500_density_0_0.8_8 -b y -cwd -V 'python test_poisson_inputs_modulate_density.py 20 2900 500 0 0.8 8'
qsub -N job_022016_e2900_i500_density_0_0.8_9 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs_modulate_density.py 20 2900 500 0 0.8 9'
