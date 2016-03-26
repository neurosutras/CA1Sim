#!/bin/bash
#cd $HOME/PycharmProjects/NEURON/
cd $HOME/CA1Sim
qsub -N job_032616_e3000_i500_comp_phase_0_inh0_0 -b y -cwd -V 'python simulate_place_cell_compressed_phase.py 0 3000 500 0 2.5 0'
qsub -N job_032616_e3000_i500_comp_phase_0_inh0_1 -b y -cwd -V 'python simulate_place_cell_compressed_phase.py 0 3000 500 0 2.5 1'
qsub -N job_032616_e3000_i500_comp_phase_0_inh0_2 -b y -cwd -V 'python simulate_place_cell_compressed_phase.py 0 3000 500 0 2.5 2'
qsub -N job_032616_e3000_i500_comp_phase_0_inh0_3 -b y -cwd -V 'python simulate_place_cell_compressed_phase.py 0 3000 500 0 2.5 3'
qsub -N job_032616_e3000_i500_comp_phase_0_inh0_4 -b y -cwd -V 'python simulate_place_cell_compressed_phase.py 0 3000 500 0 2.5 4'
qsub -N job_032616_e3000_i500_comp_phase_0_inh0_5 -b y -cwd -V 'python simulate_place_cell_compressed_phase.py 0 3000 500 0 2.5 5'
qsub -N job_032616_e3000_i500_comp_phase_0_inh0_6 -b y -cwd -V 'python simulate_place_cell_compressed_phase.py 0 3000 500 0 2.5 6'
qsub -N job_032616_e3000_i500_comp_phase_0_inh0_7 -b y -cwd -V 'python simulate_place_cell_compressed_phase.py 0 3000 500 0 2.5 7'
qsub -N job_032616_e3000_i500_comp_phase_0_inh0_8 -b y -cwd -V 'python simulate_place_cell_compressed_phase.py 0 3000 500 0 2.5 8'
qsub -N job_032616_e3000_i500_comp_phase_0_inh0_9 -b y -cwd -V 'python simulate_place_cell_compressed_phase.py 0 3000 500 0 2.5 9'
qsub -N job_032616_e3000_i500_comp_phase_0_inh2_20 -b y -cwd -V 'python simulate_place_cell_compressed_phase.py 0 3000 500 2 2.5 20'
qsub -N job_032616_e3000_i500_comp_phase_0_inh2_21 -b y -cwd -V 'python simulate_place_cell_compressed_phase.py 0 3000 500 2 2.5 21'
qsub -N job_032616_e3000_i500_comp_phase_0_inh2_22 -b y -cwd -V 'python simulate_place_cell_compressed_phase.py 0 3000 500 2 2.5 22'
qsub -N job_032616_e3000_i500_comp_phase_0_inh2_23 -b y -cwd -V 'python simulate_place_cell_compressed_phase.py 0 3000 500 2 2.5 23'
qsub -N job_032616_e3000_i500_comp_phase_0_inh2_24 -b y -cwd -V 'python simulate_place_cell_compressed_phase.py 0 3000 500 2 2.5 24'
qsub -N job_032616_e3000_i500_comp_phase_0_inh2_25 -b y -cwd -V 'python simulate_place_cell_compressed_phase.py 0 3000 500 2 2.5 25'
qsub -N job_032616_e3000_i500_comp_phase_0_inh2_26 -b y -cwd -V 'python simulate_place_cell_compressed_phase.py 0 3000 500 2 2.5 26'
qsub -N job_032616_e3000_i500_comp_phase_0_inh2_27 -b y -cwd -V 'python simulate_place_cell_compressed_phase.py 0 3000 500 2 2.5 27'
qsub -N job_032616_e3000_i500_comp_phase_0_inh2_28 -b y -cwd -V 'python simulate_place_cell_compressed_phase.py 0 3000 500 2 2.5 28'
qsub -N job_032616_e3000_i500_comp_phase_0_inh2_29 -b y -cwd -V 'python simulate_place_cell_compressed_phase.py 0 3000 500 2 2.5 29'
qsub -N job_032616_e3200_i500_comp_phase_1_inh0_30 -b y -cwd -V 'python simulate_place_cell_compressed_phase.py 1 3200 500 0 2.5 30'
qsub -N job_032616_e3200_i500_comp_phase_1_inh0_31 -b y -cwd -V 'python simulate_place_cell_compressed_phase.py 1 3200 500 0 2.5 31'
qsub -N job_032616_e3200_i500_comp_phase_1_inh0_32 -b y -cwd -V 'python simulate_place_cell_compressed_phase.py 1 3200 500 0 2.5 32'
qsub -N job_032616_e3200_i500_comp_phase_1_inh0_33 -b y -cwd -V 'python simulate_place_cell_compressed_phase.py 1 3200 500 0 2.5 33'
qsub -N job_032616_e3200_i500_comp_phase_1_inh0_34 -b y -cwd -V 'python simulate_place_cell_compressed_phase.py 1 3200 500 0 2.5 34'
qsub -N job_032616_e3200_i500_comp_phase_1_inh0_35 -b y -cwd -V 'python simulate_place_cell_compressed_phase.py 1 3200 500 0 2.5 35'
qsub -N job_032616_e3200_i500_comp_phase_1_inh0_36 -b y -cwd -V 'python simulate_place_cell_compressed_phase.py 1 3200 500 0 2.5 36'
qsub -N job_032616_e3200_i500_comp_phase_1_inh0_37 -b y -cwd -V 'python simulate_place_cell_compressed_phase.py 1 3200 500 0 2.5 37'
qsub -N job_032616_e3200_i500_comp_phase_1_inh0_38 -b y -cwd -V 'python simulate_place_cell_compressed_phase.py 1 3200 500 0 2.5 38'
qsub -N job_032616_e3200_i500_comp_phase_1_inh0_39 -b y -cwd -V 'python simulate_place_cell_compressed_phase.py 1 3200 500 0 2.5 39'
qsub -N job_032616_e3200_i500_comp_phase_1_inh2_50 -b y -cwd -V 'python simulate_place_cell_compressed_phase.py 1 3200 500 2 2.5 50'
qsub -N job_032616_e3200_i500_comp_phase_1_inh2_51 -b y -cwd -V 'python simulate_place_cell_compressed_phase.py 1 3200 500 2 2.5 51'
qsub -N job_032616_e3200_i500_comp_phase_1_inh2_52 -b y -cwd -V 'python simulate_place_cell_compressed_phase.py 1 3200 500 2 2.5 52'
qsub -N job_032616_e3200_i500_comp_phase_1_inh2_53 -b y -cwd -V 'python simulate_place_cell_compressed_phase.py 1 3200 500 2 2.5 53'
qsub -N job_032616_e3200_i500_comp_phase_1_inh2_54 -b y -cwd -V 'python simulate_place_cell_compressed_phase.py 1 3200 500 2 2.5 54'
qsub -N job_032616_e3200_i500_comp_phase_1_inh2_55 -b y -cwd -V 'python simulate_place_cell_compressed_phase.py 1 3200 500 2 2.5 55'
qsub -N job_032616_e3200_i500_comp_phase_1_inh2_56 -b y -cwd -V 'python simulate_place_cell_compressed_phase.py 1 3200 500 2 2.5 56'
qsub -N job_032616_e3200_i500_comp_phase_1_inh2_57 -b y -cwd -V 'python simulate_place_cell_compressed_phase.py 1 3200 500 2 2.5 57'
qsub -N job_032616_e3200_i500_comp_phase_1_inh2_58 -b y -cwd -V 'python simulate_place_cell_compressed_phase.py 1 3200 500 2 2.5 58'
qsub -N job_032616_e3200_i500_comp_phase_1_inh2_59 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python simulate_place_cell_compressed_phase.py 1 3200 500 2 2.5 59'