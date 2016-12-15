#!/bin/bash
cd $HOME/CA1Sim_dev
qsub -pe batch 2 -N long_cell_121216_cell19 -l d_rt=36000 -b y -cwd -V 'python optimize_long_time_scale_plasticity_rule_121016.py 19'
qsub -pe batch 2 -N long_cell_121216_cell20 -l d_rt=36000 -b y -cwd -V 'python optimize_long_time_scale_plasticity_rule_121016.py 20'
qsub -pe batch 2 -N long_cell_121216_cell21 -l d_rt=36000 -b y -cwd -V 'python optimize_long_time_scale_plasticity_rule_121016.py 21'
qsub -pe batch 2 -N long_cell_121216_cell22 -l d_rt=36000 -b y -cwd -V 'python optimize_long_time_scale_plasticity_rule_121016.py 22'
qsub -pe batch 2 -N long_cell_121216_cell23 -l d_rt=36000 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python optimize_long_time_scale_plasticity_rule_121016.py 23'

qsub -pe batch 2 -N short_cell_121216_cell19 -l d_rt=36000 -b y -cwd -V 'python optimize_short_time_scale_plasticity_rule_121016.py 19'
qsub -pe batch 2 -N short_cell_121216_cell20 -l d_rt=36000 -b y -cwd -V 'python optimize_short_time_scale_plasticity_rule_121016.py 20'
qsub -pe batch 2 -N short_cell_121216_cell21 -l d_rt=36000 -b y -cwd -V 'python optimize_short_time_scale_plasticity_rule_121016.py 21'
qsub -pe batch 2 -N short_cell_121216_cell22 -l d_rt=36000 -b y -cwd -V 'python optimize_short_time_scale_plasticity_rule_121016.py 22'
qsub -pe batch 2 -N short_cell_121216_cell23 -l d_rt=36000 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python optimize_short_time_scale_plasticity_rule_121016.py 23'