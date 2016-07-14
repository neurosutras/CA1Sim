#!/bin/bash
#cd $HOME/PycharmProjects/NEURON/
cd $HOME/CA1Sim
qsub -N job_071316_e3200_i600_no_precess5_inh3_170 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_no_precession.py 0 3200 600 3 2.5 170'
qsub -N job_071316_e3200_i600_no_precess5_inh3_171 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_no_precession.py 0 3200 600 3 2.5 171'
qsub -N job_071316_e3200_i600_no_precess5_inh3_172 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_no_precession.py 0 3200 600 3 2.5 172'
qsub -N job_071316_e3200_i600_no_precess5_inh3_173 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_no_precession.py 0 3200 600 3 2.5 173'
qsub -N job_071316_e3200_i600_no_precess5_inh3_174 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_no_precession.py 0 3200 600 3 2.5 174'
qsub -N job_071316_e3200_i600_no_precess5_inh3_175 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_no_precession.py 0 3200 600 3 2.5 175'
qsub -N job_071316_e3200_i600_no_precess5_inh3_176 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_no_precession.py 0 3200 600 3 2.5 176'
qsub -N job_071316_e3200_i600_no_precess5_inh3_177 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_no_precession.py 0 3200 600 3 2.5 177'
qsub -N job_071316_e3200_i600_no_precess5_inh3_178 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_no_precession.py 0 3200 600 3 2.5 178'
qsub -N job_071316_e3200_i600_no_precess5_inh3_179 -l d_rt=36000 -b y -cwd -V 'python simulate_place_cell_no_precession.py 0 3200 600 3 2.5 179'
