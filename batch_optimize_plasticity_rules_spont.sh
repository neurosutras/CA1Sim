#!/bin/bash
cd $HOME/CA1Sim_dev
qsub -pe batch 2 -N long_spont_012417_cell01 -l d_rt=36000 -b y -cwd -V 'python optimize_long_time_scale_plasticity_rule_spont_011817.py 1'
qsub -pe batch 2 -N long_spont_012417_cell02 -l d_rt=36000 -b y -cwd -V 'python optimize_long_time_scale_plasticity_rule_spont_011817.py 2'
qsub -pe batch 2 -N long_spont_012417_cell03 -l d_rt=36000 -b y -cwd -V 'python optimize_long_time_scale_plasticity_rule_spont_011817.py 3'
qsub -pe batch 2 -N long_spont_012417_cell04 -l d_rt=36000 -b y -cwd -V 'python optimize_long_time_scale_plasticity_rule_spont_011817.py 4'
qsub -pe batch 2 -N long_spont_012417_cell05 -l d_rt=36000 -b y -cwd -V 'python optimize_long_time_scale_plasticity_rule_spont_011817.py 5'
qsub -pe batch 2 -N long_spont_012417_cell06 -l d_rt=36000 -b y -cwd -V 'python optimize_long_time_scale_plasticity_rule_spont_011817.py 6'
qsub -pe batch 2 -N long_spont_012417_cell07 -l d_rt=36000 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python optimize_long_time_scale_plasticity_rule_spont_011817.py 7'

qsub -pe batch 2 -N short_spont_012417_cell01 -l d_rt=36000 -b y -cwd -V 'python optimize_short_time_scale_plasticity_rule_spont_011817.py 1'
qsub -pe batch 2 -N short_spont_012417_cell02 -l d_rt=36000 -b y -cwd -V 'python optimize_short_time_scale_plasticity_rule_spont_011817.py 2'
qsub -pe batch 2 -N short_spont_012417_cell03 -l d_rt=36000 -b y -cwd -V 'python optimize_short_time_scale_plasticity_rule_spont_011817.py 3'
qsub -pe batch 2 -N short_spont_012417_cell04 -l d_rt=36000 -b y -cwd -V 'python optimize_short_time_scale_plasticity_rule_spont_011817.py 4'
qsub -pe batch 2 -N short_spont_012417_cell05 -l d_rt=36000 -b y -cwd -V 'python optimize_short_time_scale_plasticity_rule_spont_011817.py 5'
qsub -pe batch 2 -N short_spont_012417_cell06 -l d_rt=36000 -b y -cwd -V 'python optimize_short_time_scale_plasticity_rule_spont_011817.py 6'
qsub -pe batch 2 -N short_spont_012417_cell07 -l d_rt=36000 -m e -M milsteina@janelia.hhmi.org -b y -cwd -V 'python optimize_short_time_scale_plasticity_rule_spont_011817.py 7'
