#!/bin/bash
cd $HOME/CA1Sim_dev
qsub -N long_cell_112716_cell01 -l d_rt=36000 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python optimize_long_time_scale_plasticity_rule_variable_run_vel_112716.py 1'
qsub -N long_cell_112716_cell02 -l d_rt=36000 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python optimize_long_time_scale_plasticity_rule_variable_run_vel_112716.py 2'
qsub -N long_cell_112716_cell03 -l d_rt=36000 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python optimize_long_time_scale_plasticity_rule_variable_run_vel_112716.py 3'
qsub -N long_cell_112716_cell04 -l d_rt=36000 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python optimize_long_time_scale_plasticity_rule_variable_run_vel_112716.py 4'
qsub -N long_cell_112716_cell05 -l d_rt=36000 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python optimize_long_time_scale_plasticity_rule_variable_run_vel_112716.py 5'
qsub -N long_cell_112716_cell06 -l d_rt=36000 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python optimize_long_time_scale_plasticity_rule_variable_run_vel_112716.py 6'
qsub -N long_cell_112716_cell07 -l d_rt=36000 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python optimize_long_time_scale_plasticity_rule_variable_run_vel_112716.py 7'
qsub -N long_cell_112716_cell08 -l d_rt=36000 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python optimize_long_time_scale_plasticity_rule_variable_run_vel_112716.py 8'
qsub -N long_cell_112716_cell09 -l d_rt=36000 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python optimize_long_time_scale_plasticity_rule_variable_run_vel_112716.py 9'
qsub -N long_cell_112716_cell10 -l d_rt=36000 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python optimize_long_time_scale_plasticity_rule_variable_run_vel_112716.py 10'
