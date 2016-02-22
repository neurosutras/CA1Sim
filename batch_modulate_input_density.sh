#!/bin/bash
cd $HOME/CA1Sim
qsub -N job_022216_e2900_i500_density_0_1.3_0 -b y -cwd -V 'python test_poisson_inputs_modulate_density.py 20 2900 500 0 1.3 0'
qsub -N job_022216_e2900_i500_density_0_1.3_1 -b y -cwd -V 'python test_poisson_inputs_modulate_density.py 20 2900 500 0 1.3 1'
qsub -N job_022216_e2900_i500_density_0_1.2_2 -b y -cwd -V 'python test_poisson_inputs_modulate_density.py 20 2900 500 0 1.2 2'
qsub -N job_022216_e2900_i500_density_0_1.2_3 -b y -cwd -V 'python test_poisson_inputs_modulate_density.py 20 2900 500 0 1.2 3'
qsub -N job_022216_e2900_i500_density_0_1.1_4 -b y -cwd -V 'python test_poisson_inputs_modulate_density.py 20 2900 500 0 1.1 4'
qsub -N job_022216_e2900_i500_density_0_1.1_5 -b y -cwd -V 'python test_poisson_inputs_modulate_density.py 20 2900 500 0 1.1 5'
qsub -N job_022216_e2900_i500_density_0_1.05_6 -b y -cwd -V 'python test_poisson_inputs_modulate_density.py 20 2900 500 0 1.05 6'
qsub -N job_022216_e2900_i500_density_0_1.05_7 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs_modulate_density.py 20 2900 500 0 1.05 7'
