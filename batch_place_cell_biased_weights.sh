#!/bin/bash
#cd $HOME/PycharmProjects/NEURON/
cd $HOME/CA1Sim
qsub -N job_052416_e3200_i500_no_na_20_inh3_620 -b y -cwd -V 'python simulate_place_cell_Type_A_no_na.py 20 3200 500 3 2.5 620'
qsub -N job_052416_e3200_i500_no_na_20_inh3_621 -b y -cwd -V 'python simulate_place_cell_Type_A_no_na.py 20 3200 500 3 2.5 621'
qsub -N job_052416_e3200_i500_no_na_20_inh3_622 -b y -cwd -V 'python simulate_place_cell_Type_A_no_na.py 20 3200 500 3 2.5 622'
qsub -N job_052416_e3200_i500_no_na_20_inh3_623 -b y -cwd -V 'python simulate_place_cell_Type_A_no_na.py 20 3200 500 3 2.5 623'
qsub -N job_052416_e3200_i500_no_na_20_inh3_624 -b y -cwd -V 'python simulate_place_cell_Type_A_no_na.py 20 3200 500 3 2.5 624'
qsub -N job_052416_e3200_i500_no_na_20_inh3_625 -b y -cwd -V 'python simulate_place_cell_Type_A_no_na.py 20 3200 500 3 2.5 625'
qsub -N job_052416_e3200_i500_no_na_20_inh3_626 -b y -cwd -V 'python simulate_place_cell_Type_A_no_na.py 20 3200 500 3 2.5 626'
qsub -N job_052416_e3200_i500_no_na_20_inh3_627 -b y -cwd -V 'python simulate_place_cell_Type_A_no_na.py 20 3200 500 3 2.5 627'
qsub -N job_052416_e3200_i500_no_na_20_inh3_628 -b y -cwd -V 'python simulate_place_cell_Type_A_no_na.py 20 3200 500 3 2.5 628'
qsub -N job_052416_e3200_i500_no_na_20_inh3_629 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python simulate_place_cell_Type_A_no_na.py 20 3200 500 3 2.5 629'
