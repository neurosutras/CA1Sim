#!/bin/bash
cd $HOME/CA1Sim
qsub -N job_021116_e3100_i500_mem_eff2_0_30 -b y -cwd -V 'python test_poisson_inputs_memory_efficient2.py 1 3100 500 0 30'
qsub -N job_021116_e3100_i500_mem_eff2_0_31 -b y -cwd -V 'python test_poisson_inputs_memory_efficient2.py 1 3100 500 0 31'
qsub -N job_021116_e3100_i500_mem_eff2_0_32 -b y -cwd -V 'python test_poisson_inputs_memory_efficient2.py 1 3100 500 0 32'
qsub -N job_021116_e3100_i500_mem_eff2_0_33 -b y -cwd -V 'python test_poisson_inputs_memory_efficient2.py 1 3100 500 0 33'
qsub -N job_021116_e3100_i500_mem_eff2_1_40 -b y -cwd -V 'python test_poisson_inputs_memory_efficient2.py 1 3100 500 1 40'
qsub -N job_021116_e3100_i500_mem_eff2_1_41 -b y -cwd -V 'python test_poisson_inputs_memory_efficient2.py 1 3100 500 1 41'
qsub -N job_021116_e3100_i500_mem_eff2_1_42 -b y -cwd -V 'python test_poisson_inputs_memory_efficient2.py 1 3100 500 1 42'
qsub -N job_021116_e3100_i500_mem_eff2_1_43 -b y -cwd -V 'python test_poisson_inputs_memory_efficient2.py 1 3100 500 1 43'
qsub -N job_021116_e3100_i500_mem_eff2_2_50 -b y -cwd -V 'python test_poisson_inputs_memory_efficient2.py 1 3100 500 2 50'
qsub -N job_021116_e3100_i500_mem_eff2_2_51 -b y -cwd -V 'python test_poisson_inputs_memory_efficient2.py 1 3100 500 2 51'
qsub -N job_021116_e3100_i500_mem_eff2_2_52 -b y -cwd -V 'python test_poisson_inputs_memory_efficient2.py 1 3100 500 2 52'
qsub -N job_021116_e3100_i500_mem_eff2_2_53 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python test_poisson_inputs_memory_efficient2.py 1 3100 500 2 53'