#!/bin/bash
#cd $HOME/PycharmProjects/NEURON/
cd $HOME/CA1Sim
qsub -N job3_021616_e6100_i500_no_nmda3_inh0_200 -b y -cwd -V 'python test_poisson_inputs_no_nmda3.py 10 6100 500 0 200'
qsub -N job3_021616_e6100_i500_no_nmda3_inh0_201 -b y -cwd -V 'python test_poisson_inputs_no_nmda3.py 10 6100 500 0 201'
qsub -N job3_021616_e6100_i500_no_nmda3_inh0_202 -b y -cwd -V 'python test_poisson_inputs_no_nmda3.py 10 6100 500 0 202'
qsub -N job3_021616_e6100_i500_no_nmda3_inh0_203 -b y -cwd -V 'python test_poisson_inputs_no_nmda3.py 10 6100 500 0 203'
qsub -N job3_021616_e6100_i500_no_nmda3_inh1_210 -b y -cwd -V 'python test_poisson_inputs_no_nmda3.py 10 6100 500 1 210'
qsub -N job3_021616_e6100_i500_no_nmda3_inh1_211 -b y -cwd -V 'python test_poisson_inputs_no_nmda3.py 10 6100 500 1 211'
qsub -N job3_021616_e6100_i500_no_nmda3_inh1_212 -b y -cwd -V 'python test_poisson_inputs_no_nmda3.py 10 6100 500 1 212'
qsub -N job3_021616_e6100_i500_no_nmda3_inh1_213 -b y -cwd -V 'python test_poisson_inputs_no_nmda3.py 10 6100 500 1 213'
qsub -N job3_021616_e6100_i500_no_nmda3_inh2_220 -b y -cwd -V 'python test_poisson_inputs_no_nmda3.py 10 6100 500 2 220'
qsub -N job3_021616_e6100_i500_no_nmda3_inh2_221 -b y -cwd -V 'python test_poisson_inputs_no_nmda3.py 10 6100 500 2 221'
qsub -N job3_021616_e6100_i500_no_nmda3_inh2_222 -b y -cwd -V 'python test_poisson_inputs_no_nmda3.py 10 6100 500 2 222'
qsub -N job3_021616_e6100_i500_no_nmda3_inh2_223 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs_no_nmda3.py 10 6100 500 2 223'