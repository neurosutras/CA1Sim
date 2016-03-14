#!/bin/bash
#cd $HOME/PycharmProjects/NEURON/
cd $HOME/CA1Sim
qsub -N job_031416_e6200_i500_no_nmda_10_inh0_204 -b y -cwd -V 'python simulate_place_cell_no_nmda.py 10 6200 500 0 204'
qsub -N job_031416_e6200_i500_no_nmda_10_inh0_205 -b y -cwd -V 'python simulate_place_cell_no_nmda.py 10 6200 500 0 205'
qsub -N job_031416_e6200_i500_no_nmda_10_inh0_206 -b y -cwd -V 'python simulate_place_cell_no_nmda.py 10 6200 500 0 206'
qsub -N job_031416_e6200_i500_no_nmda_10_inh0_207 -b y -cwd -V 'python simulate_place_cell_no_nmda.py 10 6200 500 0 207'
qsub -N job_031416_e6200_i500_no_nmda_10_inh0_208 -b y -cwd -V 'python simulate_place_cell_no_nmda.py 10 6200 500 0 208'
qsub -N job_031416_e6200_i500_no_nmda_10_inh0_209 -b y -cwd -V 'python simulate_place_cell_no_nmda.py 10 6200 500 0 209'
qsub -N job_031416_e6200_i500_no_nmda_10_inh1_214 -b y -cwd -V 'python simulate_place_cell_no_nmda.py 10 6200 500 1 214'
qsub -N job_031416_e6200_i500_no_nmda_10_inh1_215 -b y -cwd -V 'python simulate_place_cell_no_nmda.py 10 6200 500 1 215'
qsub -N job_031416_e6200_i500_no_nmda_10_inh1_216 -b y -cwd -V 'python simulate_place_cell_no_nmda.py 10 6200 500 1 216'
qsub -N job_031416_e6200_i500_no_nmda_10_inh1_217 -b y -cwd -V 'python simulate_place_cell_no_nmda.py 10 6200 500 1 217'
qsub -N job_031416_e6200_i500_no_nmda_10_inh1_218 -b y -cwd -V 'python simulate_place_cell_no_nmda.py 10 6200 500 1 218'
qsub -N job_031416_e6200_i500_no_nmda_10_inh1_219 -b y -cwd -V 'python simulate_place_cell_no_nmda.py 10 6200 500 1 219'
qsub -N job_031416_e6200_i500_no_nmda_10_inh2_224 -b y -cwd -V 'python simulate_place_cell_no_nmda.py 10 6200 500 2 224'
qsub -N job_031416_e6200_i500_no_nmda_10_inh2_225 -b y -cwd -V 'python simulate_place_cell_no_nmda.py 10 6200 500 2 225'
qsub -N job_031416_e6200_i500_no_nmda_10_inh2_226 -b y -cwd -V 'python simulate_place_cell_no_nmda.py 10 6200 500 2 226'
qsub -N job_031416_e6200_i500_no_nmda_10_inh2_227 -b y -cwd -V 'python simulate_place_cell_no_nmda.py 10 6200 500 2 227'
qsub -N job_031416_e6200_i500_no_nmda_10_inh2_228 -b y -cwd -V 'python simulate_place_cell_no_nmda.py 10 6200 500 2 228'
qsub -N job_031416_e6200_i500_no_nmda_10_inh2_229 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python simulate_place_cell_no_nmda.py 10 6200 500 2 229'

