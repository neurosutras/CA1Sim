#!/bin/bash
#cd $HOME/PycharmProjects/NEURON/
cd $HOME/CA1Sim
qsub -N job_070116_e3200_i0_subtr0_inh0_10 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh.py 0 3200 0 0 2.5 10'
qsub -N job_070116_e3200_i0_subtr0_inh0_11 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh.py 0 3200 0 0 2.5 11'
qsub -N job_070116_e3200_i0_subtr0_inh0_12 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh.py 0 3200 0 0 2.5 12'
qsub -N job_070116_e3200_i0_subtr0_inh0_13 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh.py 0 3200 0 0 2.5 13'
qsub -N job_070116_e3200_i0_subtr0_inh0_14 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh.py 0 3200 0 0 2.5 14'
qsub -N job_070116_e3200_i0_subtr0_inh0_15 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh.py 0 3200 0 0 2.5 15'
qsub -N job_070116_e3200_i0_subtr0_inh0_16 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh.py 0 3200 0 0 2.5 16'
qsub -N job_070116_e3200_i0_subtr0_inh0_17 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh.py 0 3200 0 0 2.5 17'
qsub -N job_070116_e3200_i0_subtr0_inh0_18 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh.py 0 3200 0 0 2.5 18'
qsub -N job_070116_e3200_i0_subtr0_inh0_19 -m e -M milsteina@janelia.hhmi.org -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_subtr_inh.py 0 3200 0 0 2.5 19'
