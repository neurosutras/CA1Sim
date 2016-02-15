#!/bin/bash
#cd $HOME/PycharmProjects/NEURON/
cd $HOME/CA1Sim
qsub -N job_021516_e5800_i500_no_nmda3_inh0_200 -b y -cwd -V 'python test_poisson_inputs_no_nmda3.py 10 5800 500 0 200'
qsub -N job_021516_e5800_i500_no_nmda3_inh0_201 -b y -cwd -V 'python test_poisson_inputs_no_nmda3.py 10 5800 500 0 201'
qsub -N job_021516_e6000_i500_no_nmda3_inh0_202 -b y -cwd -V 'python test_poisson_inputs_no_nmda3.py 10 6000 500 0 202'
qsub -N job_021516_e6000_i500_no_nmda3_inh0_203 -b y -cwd -V 'python test_poisson_inputs_no_nmda3.py 10 6000 500 0 203'
qsub -N job_021516_e6200_i500_no_nmda3_inh0_204 -b y -cwd -V 'python test_poisson_inputs_no_nmda3.py 10 6200 500 0 204'
qsub -N job_021516_e6200_i500_no_nmda3_inh0_205 -b y -cwd -V 'python test_poisson_inputs_no_nmda3.py 10 6200 500 0 205'
qsub -N job_021516_e6400_i500_no_nmda3_inh0_206 -b y -cwd -V 'python test_poisson_inputs_no_nmda3.py 10 6400 500 0 206'
qsub -N job_021516_e6400_i500_no_nmda3_inh0_207 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs_no_nmda3.py 10 6400 500 0 207'
qsub -N job_021516_e5600_i500_no_nmda3_inh0_208 -b y -cwd -V 'python test_poisson_inputs_no_nmda3.py 10 5600 500 0 208'
qsub -N job_021516_e5600_i500_no_nmda3_inh0_209 -b y -cwd -V 'python test_poisson_inputs_no_nmda3.py 10 5600 500 0 209'