#!/bin/bash
#cd $HOME/PycharmProjects/NEURON/
cd $HOME/CA1Sim_dev
qsub -N job_113016_e1600_i0_runvel0_0_0 -m e -M milsteina@janelia.hhmi.org -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_load_weights_wrap_around.py 0 1600 600 0'
qsub -N job_113016_e1600_i0_runvel0_0_1 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_load_weights_wrap_around.py 0 1600 600 1'
qsub -N job_113016_e1600_i0_runvel0_0_2 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_load_weights_wrap_around.py 0 1600 600 2'
qsub -N job_113016_e1600_i0_runvel0_0_3 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_load_weights_wrap_around.py 0 1600 600 3'
qsub -N job_113016_e1700_i0_runvel0_0_0 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_load_weights_wrap_around.py 0 1700 600 0'
qsub -N job_113016_e1700_i0_runvel0_0_1 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_load_weights_wrap_around.py 0 1700 600 1'
qsub -N job_113016_e1700_i0_runvel0_0_2 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_load_weights_wrap_around.py 0 1700 600 2'
qsub -N job_113016_e1700_i0_runvel0_0_3 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_load_weights_wrap_around.py 0 1700 600 3'
qsub -N job_113016_e1500_i0_runvel0_0_0 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_load_weights_wrap_around.py 0 1500 600 0'
qsub -N job_113016_e1500_i0_runvel0_0_1 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_load_weights_wrap_around.py 0 1500 600 1'
qsub -N job_113016_e1500_i0_runvel0_0_2 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_load_weights_wrap_around.py 0 1500 600 2'
qsub -N job_113016_e1500_i0_runvel0_0_3 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_load_weights_wrap_around.py 0 1500 600 3'
