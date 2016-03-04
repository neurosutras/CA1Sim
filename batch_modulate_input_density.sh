#!/bin/bash
cd $HOME/CA1Sim
qsub -N job3_030416_e2900_i500_density_20_0_1.25_0 -b y -cwd -V 'python test_poisson_inputs_modulate_density.py 20 2900 500 0 1.25 0'
qsub -N job3_030416_e2900_i500_density_20_0_1.25_1 -b y -cwd -V 'python test_poisson_inputs_modulate_density.py 20 2900 500 0 1.25 1'
qsub -N job3_030416_e2900_i500_density_20_0_1.25_2 -b y -cwd -V 'python test_poisson_inputs_modulate_density.py 20 2900 500 0 1.25 2'
qsub -N job3_030416_e2900_i500_density_20_0_1.25_3 -b y -cwd -V 'python test_poisson_inputs_modulate_density.py 20 2900 500 0 1.25 3'
qsub -N job3_030416_e2900_i500_density_20_0_1.25_4 -b y -cwd -V 'python test_poisson_inputs_modulate_density.py 20 2900 500 0 1.25 4'
qsub -N job3_030416_e2900_i500_density_20_0_1.25_5 -b y -cwd -V 'python test_poisson_inputs_modulate_density.py 20 2900 500 0 1.25 5'
qsub -N job3_030416_e2900_i500_density_20_0_1.25_6 -b y -cwd -V 'python test_poisson_inputs_modulate_density.py 20 2900 500 0 1.25 6'
qsub -N job3_030416_e2900_i500_density_20_0_1.25_7 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs_modulate_density.py 20 2900 500 0 1.25 7'
qsub -N job3_030416_e2900_i500_density_20_0_1.25_8 -b y -cwd -V 'python test_poisson_inputs_modulate_density.py 20 2900 500 0 1.25 8'
qsub -N job3_030416_e2900_i500_density_20_0_1.25_9 -b y -cwd -V 'python test_poisson_inputs_modulate_density.py 20 2900 500 0 1.25 9'
